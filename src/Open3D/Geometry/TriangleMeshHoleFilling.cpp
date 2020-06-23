// ----------------------------------------------------------------------------
// -                        Open3D: www.open3d.org                            -
// ----------------------------------------------------------------------------
// The MIT License (MIT)
//
// Copyright (c) 2020 www.open3d.org
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
// IN THE SOFTWARE.
// ----------------------------------------------------------------------------

#include "Open3D/Geometry/HalfEdgeTriangleMesh.h"
#include "Open3D/Geometry/TriangleMesh.h"
#include "Open3D/Utility/Console.h"

#include <Eigen/Dense>

#include <cmath>
#include <iostream>
#include <limits>
#include <set>
#include <string>

namespace open3d {
namespace geometry {

class Weight {
public:
    Weight(const double cosine_of_max_dihedral_angle, const double area)
        : cosine_of_max_dihedral_angle_(cosine_of_max_dihedral_angle),
          area_(area) {}

    Weight operator+(const Weight &w) {
        return Weight(this->cosine_of_max_dihedral_angle_ +
                              w.cosine_of_max_dihedral_angle_,
                      this->area_ + w.area_);
    }

    bool operator<(const Weight &w) {
        return this->cosine_of_max_dihedral_angle_ <
                       w.cosine_of_max_dihedral_angle_ ||
               (this->cosine_of_max_dihedral_angle_ ==
                        w.cosine_of_max_dihedral_angle_ &&
                this->area_ < w.area_);
    }

private:
    double cosine_of_max_dihedral_angle_;
    double area_;
};

class HoleFilling {
public:
    HoleFilling(const TriangleMesh &mesh) {
        mesh_ = std::make_shared<TriangleMesh>();
        mesh_->vertices_ = mesh.vertices_;
        mesh_->vertex_normals_ = mesh.vertex_normals_;
        mesh_->vertex_colors_ = mesh.vertex_colors_;
        mesh_->triangles_ = mesh.triangles_;

        mesh_->RemoveDuplicatedVertices();
        mesh_->RemoveDuplicatedTriangles();
        mesh_->RemoveUnreferencedVertices();
        mesh_->RemoveDegenerateTriangles();
    }

    double ComputeCosineOfDihedralAngle(
            std::vector<Eigen::Vector3d> &first_point_set,
            std::vector<Eigen::Vector3d> &second_point_set) {
        const Eigen::Vector3d n1 =
                ((first_point_set[2] - first_point_set[1])
                         .cross(first_point_set[0] - first_point_set[1]))
                        .normalized();
        const Eigen::Vector3d n2 =
                ((second_point_set[2] - second_point_set[1])
                         .cross(second_point_set[0] - second_point_set[1]))
                        .normalized();
        return n1.dot(n2);
    }

    double ComputeMaxCosineOfDihedralAngle(int vidx0, int vidx1, int vidx2) {
        // Find the adjacent triangle indices
        std::unordered_set<int> adjacent_triangle_indices(3);
        auto edges_to_triangles = mesh_->GetEdgeToTrianglesMap();

        for (auto tnb :
             edges_to_triangles[mesh_->GetOrderedEdge(vidx0, vidx1)]) {
            adjacent_triangle_indices.insert(tnb);
        }
        for (auto tnb :
             edges_to_triangles[mesh_->GetOrderedEdge(vidx0, vidx2)]) {
            adjacent_triangle_indices.insert(tnb);
        }
        for (auto tnb :
             edges_to_triangles[mesh_->GetOrderedEdge(vidx1, vidx2)]) {
            adjacent_triangle_indices.insert(tnb);
        }

        // Make a point set which spans the plane given by the three vertices
        // vidx0, vidx1, and vidx2
        std::vector<Eigen::Vector3d> first_point_set = {
                mesh_->vertices_[vidx0], mesh_->vertices_[vidx1],
                mesh_->vertices_[vidx2]};

        // Find the cosine of the highest dihedral angle
        double cosine_of_max_dihedral_angle =
                std::numeric_limits<double>::lowest();
        for (int tidx : adjacent_triangle_indices) {
            const Eigen::Vector3i triangle = mesh_->triangles_[tidx];
            std::vector<Eigen::Vector3d> second_point_set = {
                    mesh_->vertices_[triangle[0]],
                    mesh_->vertices_[triangle[1]],
                    mesh_->vertices_[triangle[2]]};
            cosine_of_max_dihedral_angle =
                    std::max(cosine_of_max_dihedral_angle,
                             ComputeCosineOfDihedralAngle(first_point_set,
                                                          second_point_set));
        }
        return cosine_of_max_dihedral_angle;
    }

    Weight ComputeWeight(int vidx0, int vidx1, int vidx2) {
        return Weight(ComputeMaxCosineOfDihedralAngle(vidx0, vidx1, vidx2),
                      TriangleMesh::ComputeTriangleArea(
                              mesh_->vertices_[vidx0], mesh_->vertices_[vidx1],
                              mesh_->vertices_[vidx2]));
    }

    void AddMinimumWeightTriangulationToMesh(
            std::vector<int> &hole,
            std::vector<std::vector<int>> &minimum_weight_index,
            int begin,
            int end) {
        if (end - begin > 1) {
            int current = minimum_weight_index[begin][end];
            mesh_->triangles_.emplace_back(
                    Eigen::Vector3i(hole[end], hole[current], hole[begin]));

            // utility::LogDebug("Add triangle {} {} {}\n",
            // std::to_string(hole[begin]), std::to_string(hole[current]),
            // std::to_string(hole[end]));

            AddMinimumWeightTriangulationToMesh(hole, minimum_weight_index,
                                                begin, current);
            AddMinimumWeightTriangulationToMesh(hole, minimum_weight_index,
                                                current, end);
        }
    }

    void IdentifyHoles() {
        // Assume that the mesh is oriented, manifold, and connected
        auto het_mesh = HalfEdgeTriangleMesh::CreateFromTriangleMesh(*mesh_);

        if (het_mesh->IsEmpty()) {
            utility::LogWarning("The mesh is non-manifold\n");
            return;
        }

        boundaries_ = het_mesh->GetBoundaries();
        utility::LogDebug("Number of holes: {}\n",
                          std::to_string(boundaries_.size()));
    }

    void TriangulateHoles() {
        for (std::vector<int> hole : boundaries_) {
            int size = int(hole.size());
            std::vector<std::vector<Weight>> minimum_weight(
                    size, std::vector<Weight>(size, Weight(0, 0)));
            std::vector<std::vector<int>> minimum_weight_index(
                    size, std::vector<int>(size, -1));

            for (size_t j = 2; j < size; j++) {
                for (size_t i = 0; i < size - j; i++) {
                    Weight weight_min =
                            Weight(std::numeric_limits<double>::max(),
                                   std::numeric_limits<double>::max());
                    int weight_index = -1;
                    int k = i + j;

                    for (size_t m = i + 1; m < k; m++) {
                        Weight weight = minimum_weight[i][m] +
                                        minimum_weight[m][k] +
                                        ComputeWeight(i, m, k);

                        if (weight < weight_min) {
                            weight_min = weight;
                            weight_index = m;
                        }
                    }

                    minimum_weight[j][k] = weight_min;
                    minimum_weight_index[i][k] = weight_index;
                }
            }

            AddMinimumWeightTriangulationToMesh(hole, minimum_weight_index, 0,
                                                size - 1);
        }
    }

    void RefineMesh() { return; }

    void FairMesh() { return; }

    std::shared_ptr<TriangleMesh> Run() {
        IdentifyHoles();
        TriangulateHoles();
        // RefineMesh();
        // FairMesh();
        return mesh_;
    }

private:
    std::shared_ptr<TriangleMesh> mesh_;
    std::vector<std::vector<int>> boundaries_;
};

std::shared_ptr<TriangleMesh> TriangleMesh::FillHoles() {
    HoleFilling hf(*this);
    return hf.Run();
};

}  // namespace geometry
}  // namespace open3d
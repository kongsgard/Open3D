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

#include "open3d/geometry/TriangleMesh.h"

#include <Eigen/Dense>
#include <cmath>
#include <iostream>
#include <limits>
#include <set>
#include <string>

#include "open3d/geometry/HalfEdgeTriangleMesh.h"
#include "open3d/utility/Console.h"

namespace open3d {
namespace geometry {

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

        edges_to_triangles_ = mesh_->GetEdgeToTrianglesMap();
    }

    std::pair<int, int> GetAdjacentHoleIndices(int index,
                                               std::vector<int> hole) {
        int hidx1 = index + 1 >= hole.size() ? 0 : index + 1;
        int hidx2 = index - 1 < 0 ? hole.size() - 1 : index - 1;
        return std::pair<int, int>(hidx1, hidx2);
    }

    double ComputeDihedralAngle(const Eigen::Vector4d &plane0,
                                const Eigen::Vector4d &plane1) {
        double dihedral_angle = std::acos(
                (plane0[0] * plane1[0] + plane0[1] * plane1[1] +
                 plane0[2] * plane1[2]) /
                (std::sqrt(plane0[0] * plane0[0] + plane0[1] * plane0[1] +
                           plane0[2] * plane0[2]) *
                 std::sqrt(plane1[0] * plane1[0] + plane1[1] * plane1[1] +
                           plane1[2] * plane1[2])));

        if (isnan(dihedral_angle)) {
            // The planes are parallel.
            return 0;
        }

        return dihedral_angle;
    }

    double GetMaxDihedralAngle(int vidx0, int vidx1, int vidx2) {
        const Eigen::Vector4d plane0 = mesh_->ComputeTrianglePlane(
                mesh_->vertices_[vidx0], mesh_->vertices_[vidx1],
                mesh_->vertices_[vidx2]);
        double max_dihedral_angle = 0;

        // Find the adjacent triangle indices
        std::unordered_set<int> adjacent_triangle_indices(2);
        for (int tnb :
             edges_to_triangles_[mesh_->GetOrderedEdge(vidx0, vidx1)]) {
            adjacent_triangle_indices.insert(tnb);
        }
        for (int tnb :
             edges_to_triangles_[mesh_->GetOrderedEdge(vidx1, vidx2)]) {
            adjacent_triangle_indices.insert(tnb);
        }
        for (int tidx : adjacent_triangle_indices) {
            Eigen::Vector4d plane1 = mesh_->GetTrianglePlane(tidx);
            double dihedral_angle = ComputeDihedralAngle(plane0, plane1);
            std::cout << "dihedral_angle: " << dihedral_angle << std::endl;
            max_dihedral_angle = std::max(max_dihedral_angle,
                                          ComputeDihedralAngle(plane0, plane1));
        }

        return max_dihedral_angle;
    }

    void IdentifyHoles() {
        // Assume that the mesh is oriented, manifold, and connected.
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
            std::vector<double> max_dihedral_angles(size, 0);
            int vidx0, vidx1, vidx2;
            for (int i = 0; i < size; i++) {
                auto adjacent_hole_indices = GetAdjacentHoleIndices(i, hole);
                max_dihedral_angles[i] = GetMaxDihedralAngle(
                        hole[adjacent_hole_indices.first], hole[i],
                        hole[adjacent_hole_indices.second]);
            }

            int hole_index;
            while (hole.size() > 3) {
                // Add a triangle for the vertex with the lowest
                // max_dihedral_angle and its two adjacent vertices in the hole.
                hole_index = std::min_element(max_dihedral_angles.begin(),
                                              max_dihedral_angles.end()) -
                             max_dihedral_angles.begin();
                auto adjacent_hole_indices =
                        GetAdjacentHoleIndices(hole_index, hole);
                mesh_->triangles_.emplace_back(Eigen::Vector3i(
                        hole[adjacent_hole_indices.first], hole[hole_index],
                        hole[adjacent_hole_indices.second]));

                // Remove the vertex with min_dihedral_angle from the hole.
                hole.erase(hole.begin() + hole_index);
                max_dihedral_angles.erase(max_dihedral_angles.begin() +
                                          hole_index);

                // Update the edges_to_triangles map.
                auto AddEdge = [&](int vidx0, int vidx1, int tidx) {
                    edges_to_triangles_[mesh_->GetOrderedEdge(vidx0, vidx1)]
                            .push_back(tidx);
                };
                int tidx = mesh_->triangles_.size() - 1;
                const auto &triangle = mesh_->triangles_[tidx];
                AddEdge(triangle(0), triangle(1), int(tidx));
                AddEdge(triangle(1), triangle(2), int(tidx));
                AddEdge(triangle(2), triangle(0), int(tidx));

                // Update minimum dihedral angle for the adjacent vertices
                adjacent_hole_indices =
                        GetAdjacentHoleIndices(hole_index, hole);
                auto adjacent_hole_indices_to_neighbor = GetAdjacentHoleIndices(
                        adjacent_hole_indices.first, hole);
                max_dihedral_angles[adjacent_hole_indices.first] =
                        GetMaxDihedralAngle(
                                adjacent_hole_indices_to_neighbor.first,
                                adjacent_hole_indices.first,
                                adjacent_hole_indices_to_neighbor.second);
                adjacent_hole_indices_to_neighbor = GetAdjacentHoleIndices(
                        adjacent_hole_indices.second, hole);
                max_dihedral_angles[adjacent_hole_indices.second] =
                        GetMaxDihedralAngle(
                                adjacent_hole_indices_to_neighbor.first,
                                adjacent_hole_indices.second,
                                adjacent_hole_indices_to_neighbor.second);
            }

            // Add the final triangle to close the hole.
            mesh_->triangles_.emplace_back(
                    Eigen::Vector3i(hole[2], hole[1], hole[0]));
        }
    }

    std::shared_ptr<TriangleMesh> Run() {
        IdentifyHoles();
        TriangulateHoles();
        return mesh_;
    }

private:
    std::shared_ptr<TriangleMesh> mesh_;
    std::vector<std::vector<int>> boundaries_;
    std::unordered_map<Eigen::Vector2i,
                       std::vector<int>,
                       utility::hash_eigen::hash<Eigen::Vector2i>>
            edges_to_triangles_;
};

std::shared_ptr<TriangleMesh> TriangleMesh::FillHoles() {
    HoleFilling hf(*this);
    return hf.Run();
};

}  // namespace geometry
}  // namespace open3d
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
#include <string>

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

        // Starting index for the triangles to be added in the triangulation
        // step
        first_unrefined_triangle_idx_ = mesh_->triangles_.size();
    }

    void IdentifyHoles() {
        // Assume that the mesh is oriented, manifold, and connected
        auto het_mesh = HalfEdgeTriangleMesh::CreateFromTriangleMesh(*mesh_);

        if (het_mesh->IsEmpty()) {
            utility::LogWarning("The mesh is non manifold\n");
            return;
        }

        boundaries_ = het_mesh->GetBoundaries();
        utility::LogDebug("Number of holes: {}\n",
                          std::to_string(boundaries_.size()));
    }

    void TriangulateHoles() {
        for (std::vector<int> hole : boundaries_) {
            int n = int(hole.size());
            std::vector<std::vector<float>> minimum_weight(
                    n, std::vector<float>(n, 0));
            std::vector<std::vector<int>> minimum_weight_index(
                    n, std::vector<int>(n, -1));

            std::vector<std::vector<std::tuple<float, float>>> minimum_weight2(
                    n, std::vector<std::tuple<float, float>>(
                               n, std::tuple<float, float>(0, 0)));

            for (size_t j = 2; j < n; j++) {
                for (size_t i = 0; i < n - j; i++) {
                    float weight_min = std::numeric_limits<float>::max();
                    int weight_index = -1;
                    int k = i + j;

                    for (size_t m = i + 1; m < k; m++) {
                        float weight = minimum_weight[i][m] +
                                       minimum_weight[m][k] +
                                       TriangleMesh::ComputeTriangleArea(
                                               mesh_->vertices_[i],
                                               mesh_->vertices_[m],
                                               mesh_->vertices_[k]);
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
                                                n - 1);
        }
    }

    double CalculateCosineOfDihedralAngle(const Eigen::Vector3d &p0,
                                          const Eigen::Vector3d &p1,
                                          const Eigen::Vector3d &p2) {
        const Eigen::Vector3d n1 = ((p2 - p1).cross(p0 - p1)).normalized();
        const Eigen::Vector3d n2 = ((p2 - p1).cross(p0 - p1)).normalized();
        return n1.dot(n2);
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

    void RefineMesh() {
        double density_control_factor = std::sqrt(2.0);
        mesh_->ComputeAdjacencyList();

        for (std::vector<int> hole : boundaries_) {
            // int n = int(hole.size());

            for (size_t tidx = first_unrefined_triangle_idx_;
                 tidx < mesh_->triangles_.size(); ++tidx) {
                const Eigen::Vector3i &triangle = mesh_->triangles_[tidx];
                const Eigen::Vector3d &vertex0 = mesh_->vertices_[triangle(0)];
                const Eigen::Vector3d &vertex1 = mesh_->vertices_[triangle(1)];
                const Eigen::Vector3d &vertex2 = mesh_->vertices_[triangle(2)];

                Eigen::Vector3d centroid(vertex0.sum() / 3, vertex1.sum() / 3,
                                         +vertex2.sum() / 3);

                // Compute the scale attribute as the average length of the
                // edges that are adjacent to the vertices in the triangle
                Eigen::Vector3d scale_attribute_vector =
                        Eigen::Vector3d::Zero();
                for (size_t vidx = 0; vidx < 3; ++vidx) {
                    for (int adjacent_vidx :
                         mesh_->adjacency_list_[triangle(vidx)]) {
                        const double edge_length =
                                ComputeVertexToVertexDistance(
                                        mesh_->vertices_[triangle(vidx)],
                                        mesh_->vertices_[adjacent_vidx]);
                        scale_attribute_vector[vidx] += edge_length;
                    }
                    scale_attribute_vector[vidx] /=
                            mesh_->adjacency_list_[triangle(vidx)].size();
                }
                double scale_attribute = scale_attribute_vector.mean();

                Eigen::Vector3d scaled_distances = Eigen::Vector3d::Zero();
                for (size_t vidx = 0; vidx < 3; ++vidx) {
                    scaled_distances[vidx] =
                            density_control_factor *
                            ComputeVertexToVertexDistance(
                                    centroid, mesh_->vertices_[triangle(vidx)]);
                }

                if ((scaled_distances[0] > scale_attribute_vector[0] &&
                     scaled_distances[0] > scale_attribute) ||
                    (scaled_distances[1] > scale_attribute_vector[1] &&
                     scaled_distances[1] > scale_attribute) ||
                    (scaled_distances[2] > scale_attribute_vector[2] &&
                     scaled_distances[2] > scale_attribute)) {
                    // Add vertex
                    int idx_centroid = mesh_->vertices_.size();
                    mesh_->vertices_.emplace_back(centroid);

                    // Remove existing triangle
                    mesh_->triangles_.erase(mesh_->triangles_.begin() + tidx);

                    // Add three new triangles
                    mesh_->triangles_.emplace_back(Eigen::Vector3i(
                            triangle(0), idx_centroid, triangle(2)));
                    mesh_->triangles_.emplace_back(Eigen::Vector3i(
                            triangle(1), idx_centroid, triangle(2)));
                    mesh_->triangles_.emplace_back(Eigen::Vector3i(
                            triangle(0), idx_centroid, triangle(1)));
                }
            }

            break;  // TODO: Remove

            // TODO: Find another way to reference the idx. This way won't work,
            // as the triangles are changed in the previous steps.
            // first_unrefined_triangle_idx_ += n;
        }
    }

    double ComputeVertexToVertexDistance(Eigen::Vector3d &vertex0,
                                         Eigen::Vector3d &vertex1) {
        return std::abs((vertex1 - vertex0).sum());
    }

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
    int first_unrefined_triangle_idx_;
};

std::shared_ptr<TriangleMesh> TriangleMesh::FillHoles() {
    HoleFilling hf(*this);
    return hf.Run();
};

}  // namespace geometry
}  // namespace open3d
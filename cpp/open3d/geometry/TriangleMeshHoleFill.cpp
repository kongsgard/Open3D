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
        mesh_->ComputeVertexNormals(/*normalized=*/false);
    }

    std::pair<int, int> GetAdjacentFrontIndices(int index,
                                                std::vector<int> hole) {
        const int hidx1 = index - 1 < 0 ? hole.size() - 1 : index - 1;
        const int hidx2 = index + 1 >= hole.size() ? 0 : index + 1;
        return std::pair<int, int>(hidx1, hidx2);
    }

    double ComputeAngle(int vidx0, int vidx1, int vidx2) {
        const Eigen::Vector3d v1 =
                mesh_->vertices_[vidx0] - mesh_->vertices_[vidx1];
        const Eigen::Vector3d v2 =
                mesh_->vertices_[vidx2] - mesh_->vertices_[vidx1];
        return std::acos(v1.dot(v2) / (v1.norm() * v2.norm()));
    }

    void ComputeAnglesForNeighbors(std::vector<double> &angles,
                                   std::vector<int> &front,
                                   std::pair<int, int> adjacent_front_indices) {
        // Update the vertex angles for vidx_prev and vidx_next
        const size_t vidx_prev = front[adjacent_front_indices.first];
        const size_t vidx_next = front[adjacent_front_indices.second];
        auto adjacent_front_indices_to_neighbors =
                GetAdjacentFrontIndices(adjacent_front_indices.first, front);
        angles[vidx_prev] =
                ComputeAngle(front[adjacent_front_indices.first], vidx_prev,
                             front[adjacent_front_indices.second]);
        adjacent_front_indices_to_neighbors =
                GetAdjacentFrontIndices(adjacent_front_indices.second, front);
        angles[vidx_prev] =
                ComputeAngle(front[adjacent_front_indices.first], vidx_next,
                             front[adjacent_front_indices.second]);
    }

    int SelectBoundaryVertexIndexByMinRandomRule(int min_angle_index,
                                                 double angle,
                                                 int num_angles) {
        if (angle < lower_angle_threshold_) {
            return min_angle_index;
        }

        // If lower_angle_threshold_ < angle < pi,
        // randomly selecting a vertex gives a better end result.
        return rand() % num_angles;
    }

    Eigen::Vector3d EstimateVertexNormal(size_t vidx,
                                         double min_angle,
                                         const Eigen::Vector3d &vertex,
                                         const Eigen::Vector3d &vertex_prev,
                                         const Eigen::Vector3d &vertex_next) {
        Eigen::Vector3d vertex_normal = mesh_->vertex_normals_[vidx];
        size_t num_adjacent_vertices = mesh_->adjacency_list_[vidx].size();

        Eigen::Vector3d v1 = vertex_prev - vertex;
        Eigen::Vector3d v2 = vertex_next - vertex;
        Eigen::Vector3d new_triangle_normal = v1.cross(v2);
        if (min_angle < lower_angle_threshold_) {
            // Do nothing
        } else if (lower_angle_threshold_ < min_angle &&
                   min_angle <= upper_angle_threshold_) {
            new_triangle_normal *= lower_angle_threshold_ / min_angle;
        } else {
            new_triangle_normal *= 2;
        }

        vertex_normal = ((num_adjacent_vertices - 1) * vertex_normal +
                         new_triangle_normal) /
                        num_adjacent_vertices;
        return vertex_normal;
    }

    bool IsConvex(size_t vidx,
                  const Eigen::Vector3d &vertex,
                  Eigen::Vector3d &vertex_normal) {
        size_t num_adjacent_vertices = mesh_->adjacency_list_[vidx].size();

        Eigen::Vector3d laplacian_coordinate = -vertex;
        // Uniform weights are used here.
        for (auto nbidx : mesh_->adjacency_list_[vidx]) {
            laplacian_coordinate +=
                    mesh_->vertices_[nbidx] / num_adjacent_vertices;
        }

        if (vertex_normal.dot(laplacian_coordinate) <= 0) {
            return true;
        }
        return false;
    }

    std::vector<Eigen::Vector3d> ComputeVertexCandidates(
            size_t min_angle_index,
            double min_angle,
            size_t vidx,
            size_t vidx_prev,
            size_t vidx_next,
            bool insert_two_vertices = false) {
        const Eigen::Vector3d &vertex = mesh_->vertices_[vidx];
        const Eigen::Vector3d &vertex_prev = mesh_->vertices_[vidx_prev];
        const Eigen::Vector3d &vertex_next = mesh_->vertices_[vidx_next];

        // Contains a bisector of the triangle [vertex_prev, vertex,
        // vertex_next] if only one vertex should be inserted.
        // If two vertices should be inserted, it contains two
        // trisectors.
        std::vector<Eigen::Vector3d> sectors;
        if (insert_two_vertices) {
            sectors.push_back((1 / 3 * (vertex_next - vertex_prev) - vertex)
                                      .normalized());
            sectors.push_back((2 / 3 * (vertex_next - vertex_prev) - vertex)
                                      .normalized());
        } else {
            sectors.push_back(
                    (0.5 * (vertex_next - vertex_prev) - vertex).normalized());
        }
        Eigen::Vector3d vertex_normal = EstimateVertexNormal(
                vidx, min_angle, vertex, vertex_prev, vertex_next);

        std::vector<Eigen::Vector3d> vertex_candidates;
        for (const auto sector : sectors) {
            const Eigen::Vector3d sector_rotated =
                    (sector - sector.dot(vertex_normal) * vertex_normal)
                            .normalized();
            const Eigen::Vector3d vertex_normal_rotated =
                    (vertex_normal +
                     alpha_ * std::abs(vertex_normal.dot(sector)) * sector)
                            .normalized();

            const double theta =
                    std::acos(vertex_normal.dot(vertex_normal_rotated));
            const double k = (std::cos(theta) - 1) /
                             vertex_normal_rotated.dot(sector_rotated);

            Eigen::Vector3d insertion_direction;
            if (IsConvex(vidx, vertex, vertex_normal)) {
                Eigen::Vector3d insertion_direction =
                        (sector_rotated + k * vertex_normal_rotated)
                                .normalized();
            } else {
                Eigen::Vector3d insertion_direction =
                        (sector_rotated - k * vertex_normal_rotated)
                                .normalized();
            }

            const Eigen::Vector3d vertex_new =
                    vertex +
                    0.5 * (vertex_prev + vertex_next - 2 * vertex).norm() *
                            insertion_direction;
            vertex_candidates.push_back(vertex_new);

            // TODO: Only add vertex_new to vertex_candidates if the smallest
            // distance between vertex_new and every vertex of front is bigger
            // than a given threshold.
        }

        return vertex_candidates;
    }

    void AddTriangle(size_t tidx0, size_t tidx1, size_t tidx2) {
        mesh_->triangles_.emplace_back(Eigen::Vector3i(tidx0, tidx1, tidx2));

        auto &triangle = mesh_->triangles_.back();
        Eigen::Vector3d v1 =
                mesh_->vertices_[triangle(1)] - mesh_->vertices_[triangle(0)];
        Eigen::Vector3d v2 =
                mesh_->vertices_[triangle(2)] - mesh_->vertices_[triangle(0)];
        mesh_->triangle_normals_.push_back(v1.cross(v2));

        mesh_->vertex_normals_[triangle(0)] += mesh_->triangle_normals_.back();
        mesh_->vertex_normals_[triangle(1)] += mesh_->triangle_normals_.back();
        mesh_->vertex_normals_[triangle(2)] += mesh_->triangle_normals_.back();
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

    // Triangulation by the modified advancing front method.
    void TriangulateHoles() {
        for (std::vector<int> front : boundaries_) {
            int size = int(front.size());
            std::vector<double> angles(size, 0);
            for (int i = 0; i < size; i++) {
                auto adjacent_front_indices = GetAdjacentFrontIndices(i, front);
                angles[i] = ComputeAngle(front[adjacent_front_indices.first],
                                         front[i],
                                         front[adjacent_front_indices.second]);
            }

            while (front.size() >= 3) {
                int min_angle_index =
                        std::min_element(angles.begin(), angles.end()) -
                        angles.begin();
                double min_angle = angles[min_angle_index];

                min_angle_index = SelectBoundaryVertexIndexByMinRandomRule(
                        min_angle_index, min_angle, angles.size());
                min_angle = angles[min_angle_index];
                const size_t vidx = front[min_angle_index];

                auto adjacent_front_indices =
                        GetAdjacentFrontIndices(min_angle_index, front);
                const size_t &vidx_prev = front[adjacent_front_indices.first];
                const size_t &vidx_next = front[adjacent_front_indices.second];

                std::vector<Eigen::Vector3d> vertex_candidates;
                if (min_angle < lower_angle_threshold_) {
                    // Do nothing.
                } else if (lower_angle_threshold_ < min_angle &&
                           min_angle <= upper_angle_threshold_) {
                    vertex_candidates =
                            ComputeVertexCandidates(min_angle_index, min_angle,
                                                    vidx, vidx_prev, vidx_next);
                } else if (upper_angle_threshold_ < min_angle &&
                           min_angle < M_PI) {
                    vertex_candidates = ComputeVertexCandidates(
                            min_angle_index, min_angle, vidx, vidx_prev,
                            vidx_next, /*insert_two_vertices=*/true);
                } else {
                    utility::LogError("Vertex angle is not valid.");
                    return;
                }

                switch (vertex_candidates.size()) {
                    case 0: {
                        // Add one triangle, but no new vertices.
                        AddTriangle(vidx_prev, vidx_next, vidx);

                        mesh_->adjacency_list_[vidx_prev].insert(vidx_next);
                        mesh_->adjacency_list_[vidx_next].insert(vidx_prev);

                        ComputeAnglesForNeighbors(angles, front,
                                                  adjacent_front_indices);
                        front.erase(front.begin() + min_angle_index);

                        break;
                    }
                    case 1: {
                        // Insert one vertex and add two triangles.
                        mesh_->vertices_.push_back(vertex_candidates.front());
                        mesh_->adjacency_list_.push_back(
                                std::unordered_set<int>{});
                        size_t vidx_new = mesh_->vertices_.size() - 1;

                        AddTriangle(vidx_prev, vidx_new, vidx);
                        AddTriangle(vidx, vidx_new, vidx_next);

                        mesh_->adjacency_list_[vidx_prev].insert(vidx_new);
                        mesh_->adjacency_list_[vidx_new].insert(vidx_prev);
                        mesh_->adjacency_list_[vidx].insert(vidx_new);
                        mesh_->adjacency_list_[vidx_new].insert(vidx);
                        mesh_->adjacency_list_[vidx_next].insert(vidx_new);
                        mesh_->adjacency_list_[vidx_new].insert(vidx_next);

                        ComputeAnglesForNeighbors(angles, front,
                                                  adjacent_front_indices);
                        front[min_angle_index] = vidx_new;
                        angles[min_angle_index] =
                                ComputeAngle(vidx_prev, vidx_new, vidx_next);

                        break;
                    }
                    case 2: {
                        // Insert two vertices and add three triangles.
                        mesh_->vertices_.push_back(vertex_candidates.front());
                        mesh_->vertices_.push_back(vertex_candidates.back());
                        mesh_->adjacency_list_.push_back(
                                std::unordered_set<int>{});
                        mesh_->adjacency_list_.push_back(
                                std::unordered_set<int>{});
                        std::pair<size_t, size_t> vidxs_new(
                                mesh_->vertices_.size() - 2,
                                mesh_->vertices_.size() - 1);

                        AddTriangle(vidx_prev, vidxs_new.first, vidx);
                        AddTriangle(vidxs_new.first, vidxs_new.second, vidx);
                        AddTriangle(vidxs_new.second, vidx_next, vidx);

                        mesh_->adjacency_list_[vidx_prev].insert(
                                vidxs_new.first);
                        mesh_->adjacency_list_[vidxs_new.first].insert(
                                vidx_prev);
                        mesh_->adjacency_list_[vidx].insert(vidxs_new.first);
                        mesh_->adjacency_list_[vidxs_new.first].insert(vidx);
                        mesh_->adjacency_list_[vidxs_new.first].insert(
                                vidxs_new.second);
                        mesh_->adjacency_list_[vidxs_new.second].insert(
                                vidxs_new.first);
                        mesh_->adjacency_list_[vidx].insert(vidxs_new.second);
                        mesh_->adjacency_list_[vidxs_new.second].insert(vidx);
                        mesh_->adjacency_list_[vidx_next].insert(
                                vidxs_new.second);
                        mesh_->adjacency_list_[vidxs_new.second].insert(
                                vidx_next);

                        ComputeAnglesForNeighbors(angles, front,
                                                  adjacent_front_indices);
                        front[min_angle_index] = vidxs_new.first;
                        angles[min_angle_index] = ComputeAngle(
                                vidx_prev, vidxs_new.first, vidxs_new.second);
                        front.insert(front.begin() + min_angle_index + 1,
                                     vidxs_new.second);
                        angles.insert(
                                angles.begin() + min_angle_index + 1,
                                ComputeAngle(vidxs_new.first, vidxs_new.second,
                                             vidx_next));
                        // TODO: Insertion into a vector is inefficient.
                        // Consider using e.g. a list.

                        break;
                    }
                    default: {
                        utility::LogError(
                                "Computed more than two vertex candidates for "
                                "triangulation.");
                        return;
                    }
                }
            }
        }
    }

    std::shared_ptr<TriangleMesh> Run() {
        IdentifyHoles();
        TriangulateHoles();

        mesh_->ComputeVertexNormals();
        return mesh_;
    }

private:
    std::shared_ptr<TriangleMesh> mesh_;
    std::vector<std::vector<int>> boundaries_;

    double lower_angle_threshold_ = M_PI * 85 / 180;
    double upper_angle_threshold_ = M_PI * 135 / 180;
    double alpha_ = 0.45;  // Can be in range [0,1]
};

std::shared_ptr<TriangleMesh> TriangleMesh::FillHoles() {
    HoleFilling hf(*this);
    return hf.Run();
};

}  // namespace geometry
}  // namespace open3d
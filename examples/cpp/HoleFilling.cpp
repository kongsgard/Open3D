#include <Eigen/Dense>
#include <iostream>

#include "open3d/Open3D.h"

int main(int argc, char *argv[]) {
    using namespace open3d;
    utility::SetVerbosityLevel(utility::VerbosityLevel::Debug);

    //*
    // auto mesh = io::CreateMeshFromFile(
    //        "/Users/sondrebk/Documents/Open3d/examples/TestData/Bunny.ply");
    auto mesh = io::CreateMeshFromFile(
            "/Users/sondrebk/Desktop/archive/Badger.off");
    // io::WriteTriangleMesh("/Users/sondrebk/Desktop/archive/Badger.off",
    // *mesh,
    //                       false, false, false, false);
    utility::LogDebug("Triangles: {}\n",
                      std::to_string(mesh->triangles_.size()));
    auto filled_mesh = mesh->FillHoles();
    //*/

    /*
    geometry::TriangleMesh mesh;
    mesh.vertices_ = {{0, 0, 0},      {1, 0, 0},   {1, 1, 0},  {0, 1, 0},
                      {0.5, 0.5, -1}, {0.5, 0, 1}, {0.5, 1, 1}};
    mesh.triangles_ = {{0, 4, 1}, {1, 4, 2}, {2, 4, 3},
                       {3, 4, 0}, {0, 1, 5}, {2, 3, 6}};
    auto filled_mesh = mesh.FillHoles();
    //*/

    /*
    geometry::TriangleMesh mesh;
    mesh.vertices_ = {{0, 0, 0},      {1, 0, 0},   {1, 1, 0},   {0, 1, 0},
                      {0.5, 0.5, -1}, {0.5, 0, 1}, {0.5, 1, 1}, {0.5, 0.5, 2}};
    mesh.triangles_ = {{0, 4, 1}, {1, 4, 2}, {2, 4, 3}, {3, 4, 0},
                       {0, 1, 5}, {2, 3, 6}, {1, 2, 7}, {3, 0, 5}};
    auto filled_mesh = mesh.FillHoles();
    //*/

    utility::LogDebug("Triangles: {}\n",
                      std::to_string(filled_mesh->triangles_.size()));

    utility::LogDebug("IsEdgeManifold: {}\n",
                      std::to_string(filled_mesh->IsEdgeManifold(false)));
    utility::LogDebug("IsVertexManifold: {}\n",
                      std::to_string(filled_mesh->IsVertexManifold()));
    utility::LogDebug("IsWatertight: {}\n",
                      std::to_string(filled_mesh->IsWatertight()));

    std::vector<std::shared_ptr<const geometry::Geometry>> meshes;
    meshes.push_back(mesh);
    // meshes.push_back(filled_mesh);
    visualization::DrawGeometries(meshes);
}
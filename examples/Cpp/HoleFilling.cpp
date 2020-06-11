#include <Eigen/Dense>
#include <iostream>

#include "Open3D/Open3D.h"

int main(int argc, char *argv[]) {
    using namespace open3d;
    utility::SetVerbosityLevel(utility::VerbosityLevel::Debug);

    auto mesh = io::CreateMeshFromFile(
            "/Users/sondrebk/Documents/Open3d/examples/TestData/Bunny.ply");

    auto filled_mesh = mesh->FillHoles();

    utility::LogDebug("IsEdgeManifold: {}\n",
                      std::to_string(filled_mesh->IsEdgeManifold(false)));
    utility::LogDebug("IsVertexManifold: {}\n",
                      std::to_string(filled_mesh->IsVertexManifold()));
    utility::LogDebug("IsWatertight: {}\n",
                      std::to_string(filled_mesh->IsWatertight()));

    std::vector<std::shared_ptr<const geometry::Geometry>> meshes;
    meshes.push_back(filled_mesh);
    visualization::DrawGeometries(meshes);
}

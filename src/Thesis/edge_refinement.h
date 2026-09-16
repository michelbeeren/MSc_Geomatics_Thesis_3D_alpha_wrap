//
// Created by Michel Beeren on 26/01/2026.
//

#ifndef THESIS_EDGE_REFINEMENT_H
#define THESIS_EDGE_REFINEMENT_H

Mesh offset_mesh(const Mesh& mesh, const double offset);
Mesh offset_mesh_by_5cm(Mesh& mesh, const double offset_);
Mesh add_midpoint_distance_tag(const Mesh& wrapped_mesh, const Mesh& original_mesh, const double filter);
Mesh refine_round_edges(const Mesh& tagged_mesh, const Mesh& original_mesh);
Mesh refine_round_edges3(const Mesh& tagged_mesh, const Mesh& original_mesh);
Mesh refine_round_edges2(const Mesh& tagged_mesh);

#endif //THESIS_EDGE_REFINEMENT_H
#pragma once

#include<iostream>
#include<queue>
#include<stack>
#include<cstring>
#include<string>
#include<complex> // for using complex numbers
#include<cfloat>

#include "zrlib_util.h" // include util package in my lib.s
#include "zrlib_ds.h"

#define OPENMESH_SUPPORT 1

#if OPENMESH_SUPPORT
#include "util/OpenMeshUtils.h"
#include <OpenMesh\Core\Mesh\PolyMeshT.hh>
#include <OpenMesh\Core\Mesh\TriMeshT.hh>
// essential tools
void zr_mesh_transform(HEMesh& mesh, float x, float y, float z) {// simlifed version quick test 
	HEMesh::VertexIter v_it, v_end(mesh.vertices_end());
	HEMesh::Point new_point;
	
	for (v_it = mesh.vertices_begin(); v_it != v_end; ++v_it) {
		// move based on current position 
		new_point[0] = mesh.point(*v_it)[0] + x;
		new_point[1] = mesh.point(*v_it)[1] + y;
		new_point[2] = mesh.point(*v_it)[2] + z;
		mesh.set_point(*v_it, new_point);
	}
}
void zr_mesh_movecenterto(HEMesh& mesh, float x, float y, float z) {
	HEMesh::VertexIter v_it, v_end(mesh.vertices_end());
	HEMesh::Point cur_center;
	cur_center[0] = 0;
	cur_center[1] = 0;
	cur_center[2] = 0;
	int num = 0;
	for (v_it = mesh.vertices_begin(); v_it != v_end; ++v_it) {
		cur_center[0] += mesh.point(*v_it)[0];
		cur_center[1] += mesh.point(*v_it)[1];
		cur_center[2] += mesh.point(*v_it)[2];
		num++;
	}
	cur_center[0] = cur_center[0] / num;
	cur_center[1] = cur_center[1] / num;
	cur_center[2] = cur_center[2] / num;

	float move_vec_x = x - cur_center[0];
	float move_vec_y = y - cur_center[1];
	float move_vec_z = z - cur_center[2];

	zr_mesh_transform(mesh, move_vec_x, move_vec_y, move_vec_z);
} 
unsigned int zr_compute_vertex_range(HEMesh& mesh) {
	HEMesh::VertexIter v_end(mesh.vertices_end());
	int numofvertices = 0;
	for (auto v_it = mesh.vertices_begin(); v_it != v_end; ++v_it) {
		numofvertices++;
	}
	return numofvertices;
}

// creation 
void zr_create_torus(HEMesh& mesh, float r, float R, int nsides, int rings) {
	assert(nsides >= 3 && rings >= 3);

	mesh.clear();

	int n = rings * nsides;
	std::vector<OpenMesh::VertexHandle> vhandles(n);
	int k = 0;
	for (int i = 0; i < rings; i++) {
		float angle1 = (float)(i * 2.0 * 3.14159 / (rings));
		OpenMesh::Vec3f center(cos(angle1) * R, 0.0f, sin(angle1) * R);
		OpenMesh::Vec3f t1(cos(angle1), 0.0, sin(angle1));
		OpenMesh::Vec3f t2(0.0f, 1.0f, 0.0f);

		for (int j = 0; j < nsides; j++) {
			float angle2 = (float)(j * 2.0 * 3.14159 / (nsides));
			vhandles[k] = mesh.add_vertex(center + (float)(sin(angle2) * r) * t1 +
				(float)(cos(angle2) * r) * t2);
			k++;
		}
	}

	for (int i = 0; i < rings; i++) {
		for (int j = 0; j < nsides; j++) {
			int a, b, c, d;
			a = (i + 1) % (rings) * (nsides)+j;
			b = (i + 1) % (rings) * (nsides)+(j + 1) % (nsides);
			c = i * (nsides)+(j + 1) % (nsides);
			d = i * (nsides)+j;
			mesh.add_face(vhandles[a], vhandles[b], vhandles[c], vhandles[d]);
		}
	}
}
void zr_create_cube(HEMesh& mesh) {
	mesh.clear();
	/* Task 2.2.1 */
	std::vector<OpenMesh::VertexHandle> vhandles(8);

	vhandles[0] = mesh.add_vertex(OpenMesh::Vec3f(-0.5f, -0.5f, -0.5f));
	vhandles[1] = mesh.add_vertex(OpenMesh::Vec3f(0.5f, -0.5f, -0.5f));
	vhandles[2] = mesh.add_vertex(OpenMesh::Vec3f(0.5f, 0.5f, -0.5f));
	vhandles[3] = mesh.add_vertex(OpenMesh::Vec3f(-0.5f, 0.5f, -0.5f));

	vhandles[4] = mesh.add_vertex(OpenMesh::Vec3f(-0.5f, -0.5f, 0.5f));
	vhandles[5] = mesh.add_vertex(OpenMesh::Vec3f(0.5f, -0.5f, 0.5f));
	vhandles[6] = mesh.add_vertex(OpenMesh::Vec3f(0.5f, 0.5f, 0.5f));
	vhandles[7] = mesh.add_vertex(OpenMesh::Vec3f(-0.5f, 0.5f, 0.5f));

	mesh.add_face(vhandles[4], vhandles[5], vhandles[6], vhandles[7]);
	mesh.add_face(vhandles[5], vhandles[1], vhandles[2], vhandles[6]);
	mesh.add_face(vhandles[1], vhandles[0], vhandles[3], vhandles[2]);
	mesh.add_face(vhandles[0], vhandles[4], vhandles[7], vhandles[3]);
	mesh.add_face(vhandles[5], vhandles[4], vhandles[0], vhandles[1]);
	mesh.add_face(vhandles[7], vhandles[6], vhandles[2], vhandles[3]);
}
void zr_create_cube_at_posi(HEMesh& mesh, float scale, OpenMesh::Vec3f translation) {
	mesh.clear();
	/* Task 2.2.1 */
	std::vector<OpenMesh::VertexHandle> vhandles(8);

	vhandles[0] = mesh.add_vertex(scale * OpenMesh::Vec3f(-0.5f, -0.5f, -0.5f) + translation);
	vhandles[1] = mesh.add_vertex(scale * OpenMesh::Vec3f(0.5f, -0.5f, -0.5f) + translation);
	vhandles[2] = mesh.add_vertex(scale * OpenMesh::Vec3f(0.5f, 0.5f, -0.5f) + translation);
	vhandles[3] = mesh.add_vertex(scale * OpenMesh::Vec3f(-0.5f, 0.5f, -0.5f) + translation);

	vhandles[4] = mesh.add_vertex(scale * OpenMesh::Vec3f(-0.5f, -0.5f, 0.5f) + translation);
	vhandles[5] = mesh.add_vertex(scale * OpenMesh::Vec3f(0.5f, -0.5f, 0.5f) + translation);
	vhandles[6] = mesh.add_vertex(scale * OpenMesh::Vec3f(0.5f, 0.5f, 0.5f) + translation);
	vhandles[7] = mesh.add_vertex(scale * OpenMesh::Vec3f(-0.5f, 0.5f, 0.5f) + translation);

	mesh.add_face(vhandles[4], vhandles[5], vhandles[6], vhandles[7]);
	mesh.add_face(vhandles[5], vhandles[1], vhandles[2], vhandles[6]);
	mesh.add_face(vhandles[1], vhandles[0], vhandles[3], vhandles[2]);
	mesh.add_face(vhandles[0], vhandles[4], vhandles[7], vhandles[3]);
	mesh.add_face(vhandles[5], vhandles[4], vhandles[0], vhandles[1]);
	mesh.add_face(vhandles[7], vhandles[6], vhandles[2], vhandles[3]);
}
void zr_create_tetrahedron(HEMesh& mesh, float a) {
	mesh.clear();

	std::vector<OpenMesh::VertexHandle> vhandles(4);

	vhandles[0] = mesh.add_vertex(OpenMesh::Vec3f(sqrt(3.0f) * a / 3.0f, 0, 0));
	vhandles[1] =
		mesh.add_vertex(OpenMesh::Vec3f(-sqrt(3.0f) * a / 6.0f, -a / 2.0f, 0));
	vhandles[2] =
		mesh.add_vertex(OpenMesh::Vec3f(-sqrt(3.0f) * a / 6.0f, a / 2.0f, 0));
	vhandles[3] = mesh.add_vertex(OpenMesh::Vec3f(0, 0, sqrt(6.0f) * a / 3.0f));

	mesh.add_face(vhandles[0], vhandles[1], vhandles[2]);
	mesh.add_face(vhandles[0], vhandles[2], vhandles[3]);
	mesh.add_face(vhandles[0], vhandles[3], vhandles[1]);
	mesh.add_face(vhandles[3], vhandles[2], vhandles[1]);
}
void zr_create_disk(HEMesh& mesh, float radius, int slices) {
	mesh.clear();

	std::vector<OpenMesh::VertexHandle> vhandles(slices + 1);
	vhandles[0] = mesh.add_vertex(OpenMesh::Vec3f(0, 0, 0));
	for (int i = 0; i < slices; i++) {
		float angle = -i * 2 * 3.14159f / slices;
		vhandles[i + 1] = mesh.add_vertex(
			OpenMesh::Vec3f(sin(angle) * radius, cos(angle) * radius, 0));
	}
	for (int i = 0; i < slices; i++)
		mesh.add_face(vhandles[0], vhandles[1 + i % slices],
			vhandles[1 + (1 + i) % slices]);
}
void zr_create_cylinder(HEMesh& mesh, float radius, float height, int stacks,int slices) {
	assert(slices >= 3 && stacks >= 1);

	mesh.clear();

	int n = 2 + slices * (stacks + 1);
	std::vector<OpenMesh::VertexHandle> vhandles(n);

	vhandles[0] = mesh.add_vertex(OpenMesh::Vec3f(0.0f, height, 0.0f));

	int k = 1;
	for (int i = 0; i < stacks + 1; i++) {
		float h = (stacks - i) * height / (stacks);

		for (int j = 0; j < slices; j++) {
			float angle2 = j * 2.0f * 3.14159f / (float)(slices);
			vhandles[k] = mesh.add_vertex(
				OpenMesh::Vec3f(cos(angle2) * radius, h, sin(angle2) * radius));

			k++;
		}
	}

	vhandles[k] = mesh.add_vertex(OpenMesh::Vec3f(0.0f, 0.0f, 0.0f));

	for (int i = 0; i < slices; i++) {
		mesh.add_face(vhandles[0], vhandles[1 + (1 + i) % slices],
			vhandles[1 + i % slices]);

		for (int j = 0; j < stacks; j++) {
			int a, b, c, d;
			a = 1 + j * slices + (i) % slices;
			b = 1 + j * slices + (1 + i) % slices;
			c = 1 + (j + 1) * slices + (1 + i) % slices;
			d = 1 + (j + 1) * slices + (i) % slices;
			mesh.add_face(vhandles[a], vhandles[b], vhandles[c], vhandles[d]);
		}
		mesh.add_face(vhandles[vhandles.size() - 1],
			vhandles[1 + (stacks)*slices + (i) % slices],
			vhandles[1 + (stacks)*slices + (1 + i) % slices]);
	}
}
void zr_create_sphere(HEMesh& mesh, float radius, int slices, int stacks) {
	assert(slices >= 3 && stacks >= 3);

	mesh.clear();

	int n = slices * (stacks - 1) + 2;
	std::vector<OpenMesh::VertexHandle> vhandles(n);

	vhandles[0] = mesh.add_vertex(OpenMesh::Vec3f(0.0f, radius, 0.0f));

	int k = 1;
	for (int i = 1; i < stacks; i++) {
		float angle1 = 3.14159f / 2.0f - i * 3.14159f / (float)stacks;
		float r = cos(angle1) * radius;
		float height = sin(angle1) * radius;

		for (int j = 0; j < slices; j++) {
			float angle2 = j * 2.0f * 3.14159f / (float)(slices);
			vhandles[k] = mesh.add_vertex(
				OpenMesh::Vec3f(cos(angle2) * r, height, sin(angle2) * r));
			k++;
		}
	}

	vhandles[k] = mesh.add_vertex(OpenMesh::Vec3f(0.0f, -radius, 0.0f));

	for (int i = 0; i < slices; i++) {
		mesh.add_face(vhandles[0], vhandles[1 + (1 + i) % slices],
			vhandles[1 + i % slices]);

		for (int j = 0; j < stacks - 2; j++) {
			int a, b, c, d;
			a = 1 + j * slices + (i) % slices;
			b = 1 + j * slices + (1 + i) % slices;
			c = 1 + (j + 1) * slices + (1 + i) % slices;
			d = 1 + (j + 1) * slices + (i) % slices;
			mesh.add_face(vhandles[a], vhandles[b], vhandles[c], vhandles[d]);
		}
		mesh.add_face(vhandles[1 + slices * (stacks - 1)],
			vhandles[1 + (stacks - 2) * slices + (i) % slices],
			vhandles[1 + (stacks - 2) * slices + (1 + i) % slices]);
	}
}
void zr_create_icosahedron(HEMesh& mesh, float radius) {
	mesh.clear();

	float a = (float)(radius * 4.0 / sqrt(10.0 + 2.0 * sqrt(5.0)));
	float h = (float)cos(2.0 * asin(a / (2.0 * radius))) * radius;
	float r2 = (float)sqrt(radius * radius - h * h);

	std::vector<OpenMesh::VertexHandle> vhandles(12);
	int k = 0;
	vhandles[k++] = mesh.add_vertex(OpenMesh::Vec3f(0, radius, 0));

	for (int i = 0; i < 5; i++)
		vhandles[k++] = mesh.add_vertex(
			OpenMesh::Vec3f((float)cos(i * 72.0 * 3.14159 / 180.0) * r2, h,
				-(float)sin(i * 72.0 * 3.14159 / 180.0) * r2));

	for (int i = 0; i < 5; i++)
		vhandles[k++] = mesh.add_vertex(OpenMesh::Vec3f(
		(float)cos(36.0 * 3.14159 / 180.0 + i * 72.0 * 3.14159 / 180.0) * r2,
			-h,
			-(float)sin(36.0 * 3.14159 / 180.0 + i * 72.0 * 3.14159 / 180.0) * r2));

	vhandles[k] = mesh.add_vertex(OpenMesh::Vec3f(0, -radius, 0));

	for (int i = 0; i < 5; i++) {
		mesh.add_face(vhandles[0], vhandles[i + 1], vhandles[(i + 1) % 5 + 1]);
		mesh.add_face(vhandles[11], vhandles[(i + 1) % 5 + 6], vhandles[i + 6]);
		mesh.add_face(vhandles[i + 1], vhandles[i + 6], vhandles[(i + 1) % 5 + 1]);
		mesh.add_face(vhandles[(i + 1) % 5 + 1], vhandles[i + 6],
			vhandles[(i + 1) % 5 + 6]);
	}
}
void zr_create_octahedron(HEMesh& mesh, float radius) {
	// creates an octaeder mesh
	// radius is the radius of the circum sphere
	mesh.clear();

	std::vector<OpenMesh::VertexHandle> vhandles(6);
	int k = 0;
	vhandles[k++] = mesh.add_vertex(OpenMesh::Vec3f(0, radius, 0));

	for (int i = 0; i < 4; i++)
		vhandles[k++] = mesh.add_vertex(
			OpenMesh::Vec3f((float)cos(i * 3.14159 / 2.0) * radius, 0,
				-(float)sin(i * 3.14159 / 2.0) * radius));

	vhandles[k++] = mesh.add_vertex(OpenMesh::Vec3f(0, -radius, 0));

	for (int i = 0; i < 4; i++) {
		mesh.add_face(vhandles[0], vhandles[i + 1], vhandles[(i + 1) % 4 + 1]);
		mesh.add_face(vhandles[5], vhandles[(i + 1) % 4 + 1], vhandles[i + 1]);
	}
}
void zr_create_unitarrow(HEMesh& mesh, float stem_radius, float head_radius,float stem_height, int slices, int stem_stacks) {
	assert(slices >= 3 && stem_stacks >= 1 && stem_height <= 1 &&
		stem_height >= 0);

	float head_height = (float)1 - stem_height;
	mesh.clear();

	int n = 2 + slices * (stem_stacks + 2);
	std::vector<OpenMesh::VertexHandle> vhandles(n);

	float height = stem_height + head_height;
	vhandles[0] = mesh.add_vertex(OpenMesh::Vec3f(0.0f, height, 0.0f));

	int k = 1;
	for (int j = 0; j < slices; j++) {
		float angle2 = j * 2.0f * 3.14159f / (float)(slices);
		vhandles[k] = mesh.add_vertex(OpenMesh::Vec3f(
			cos(angle2) * head_radius, stem_height, sin(angle2) * head_radius));
		k++;
	}

	for (int i = 0; i < stem_stacks + 1; i++) {
		float h = (stem_stacks - i) * stem_height / (stem_stacks);

		for (int j = 0; j < slices; j++) {
			float angle2 = j * 2.0f * 3.14159f / (float)(slices);
			vhandles[k] = mesh.add_vertex(OpenMesh::Vec3f(
				cos(angle2) * stem_radius, h, sin(angle2) * stem_radius));
			k++;
		}
	}
	vhandles[k] = mesh.add_vertex(OpenMesh::Vec3f(0.0f, 0.0f, 0.0f));

	for (int i = 0; i < slices; i++) {
		mesh.add_face(vhandles[0], vhandles[1 + (1 + i) % slices],
			vhandles[1 + i % slices]);

		for (int j = 0; j < stem_stacks + 1; j++) {
			int a, b, c, d;
			a = 1 + j * slices + (i) % slices;
			b = 1 + j * slices + (1 + i) % slices;
			c = 1 + (j + 1) * slices + (1 + i) % slices;
			d = 1 + (j + 1) * slices + (i) % slices;
			mesh.add_face(vhandles[a], vhandles[b], vhandles[c], vhandles[d]);
		}
		mesh.add_face(vhandles[vhandles.size() - 1],
			vhandles[1 + (stem_stacks + 1) * slices + (i) % slices],
			vhandles[1 + (stem_stacks + 1) * slices + (1 + i) % slices]);
	}
}

// extraction
unsigned int zr_extract_shells(HEMesh& m, OpenMesh::FPropHandleT<int> perFaceShellIndex) {
	// reset the shell indices to -1 for every face
	for (auto f : m.faces()) {
		m.property(perFaceShellIndex, f) = 0;
	};
	/*Task 2.2.3*/
	// init to component/set per vertex
	// growing from one vertex

	int numOfComps = 0;

	// init
	OpenMesh::VPropHandleT<bool> visited;
	m.add_property(visited);
	HEMesh::VertexIter v_end = m.vertices_end();
	HEMesh::VertexIter v_it;
	for (v_it = m.vertices_begin(); v_it != v_end; ++v_it)
		m.property(visited, *v_it) = false;
	HEMesh::VertexHandle vh;
	HEMesh::VertexIter current_pos = m.vertices_begin();

	while (true) {
		// find an unvisited vertex
		bool found = false;
		for (v_it = current_pos; v_it != v_end; ++v_it) {
			if (!m.property(visited, *v_it)) {
				found = true;
				vh = *v_it;
				m.property(visited, *v_it) = true;
				current_pos = v_it;
				break; // do this as soon as a unvisited vertex is found
			}
		}

		// if none was found -> finished
		if (!found)
			break;
		numOfComps++;
		// here,new com found
		std::vector<HEMesh::VertexHandle> handles;
		handles.push_back(vh);

		// grow from found vertex
		while (!handles.empty()) {
			HEMesh::VertexHandle current = handles.back();
			handles.pop_back(); // getting the last one in the list
			HEMesh::VertexVertexIter vv_it;

			for (vv_it = m.vv_iter(current); vv_it.is_valid(); ++vv_it)
				if (!m.property(visited, *vv_it)) {
					m.property(visited, *vv_it) = true;
					for (HEMesh::FaceHandle vf_it : m.vf_range(vv_it)) {
						m.property(perFaceShellIndex, vf_it) = numOfComps;
					}
					handles.push_back(*vv_it);
				}
		}
	}
	m.remove_property(visited);
	return numOfComps;
}
void zr_grow_with_given_point(HEMesh& mesh, OpenMesh::FPropHandleT<int> perFaceShellIndex, unsigned int seed_vertex_idx) {
	// reset the shell indices to -1 for every face
	for (auto f : mesh.faces()) {
		mesh.property(perFaceShellIndex, f) = 0;
	};

	OpenMesh::VPropHandleT<bool> visited;
	mesh.add_property(visited);
	HEMesh::VertexIter v_end = mesh.vertices_end();
	HEMesh::VertexIter v_it = mesh.vertices_begin();
	for (v_it = mesh.vertices_begin(); v_it != v_end; ++v_it)
		mesh.property(visited, *v_it) = false;
	HEMesh::VertexHandle vh;

	int current_pos = 0;
	for (v_it = mesh.vertices_begin(); v_it != v_end; ++v_it) {
		current_pos++;
		if (current_pos == seed_vertex_idx)
		{
			vh = *v_it;
			break;
		}
	}

	std::vector<HEMesh::VertexHandle> handles;
	handles.push_back(vh);

	// grow from found vertex
	while (!handles.empty()) {
		HEMesh::VertexHandle current = handles.back();
		handles.pop_back(); // getting the last one in the list
		HEMesh::VertexVertexIter vv_it;

		for (vv_it = mesh.vv_iter(current); vv_it.is_valid(); ++vv_it)
			if (!mesh.property(visited, *vv_it)) {
				mesh.property(visited, *vv_it) = true;
				for (HEMesh::FaceHandle vf_it : mesh.vf_range(vv_it)) {
					mesh.property(perFaceShellIndex, vf_it) = 1;
				}
				handles.push_back(*vv_it);
			}
	}
}
unsigned int zr_extract_tri_strips(HEMesh& mesh, OpenMesh::FPropHandleT<int> perFaceStripIdProperty, unsigned int nTrials) {
	// prepare random engine
	std::mt19937 eng;
	std::normal_distribution<float> dist;
	int face_num = 0;
	OpenMesh::FPropHandleT<bool> faceUsedProperty;
	mesh.add_property(faceUsedProperty);
	// initialize strip index to -1 for each face
	for (auto f : mesh.faces()) {
		mesh.property(perFaceStripIdProperty, f) = -1;
		mesh.property(faceUsedProperty, f) = false;
		++face_num;
	}

	int nStrips = 0;

	/*Task 2.2.5*/

	// 1th , sample a random tri face form mesh
	// add all faces

	sample_set<HEMesh::FaceHandle> sampleset;
	sampleset.reserve(face_num);
	for (HEMesh::FaceIter f_it = mesh.faces_begin(); f_it != mesh.faces_end();
		++f_it) {
		sampleset.insert(*f_it);
	}
	std::cout << "size now:" << sampleset.size() << std::endl;

	bool not_finish = true;
	while (not_finish) {
		// check if all faces are visited
		not_finish = false;
		for (auto f : mesh.faces()) {
			if (mesh.property(perFaceStripIdProperty, f) == -1) {
				not_finish = true;
				break;
			}
		}
		if (!not_finish)
			break;
		// sample nTrials faces each time
		for (int k = 0; k < nTrials; k++) {
			HEMesh::FaceHandle k_facehandle = sampleset.sample(eng);
			std::cout << "face handle selected:" << k_facehandle << std::endl;

			// choose the longest for each point selected
			// three edges in the face
			HEMesh::HalfedgeHandle h[3];
			h[0] = mesh.halfedge_handle(k_facehandle);
			h[1] = mesh.next_halfedge_handle(h[0]);
			h[2] = mesh.next_halfedge_handle(h[1]);
			std::cout << "halfedge handle selected:" << h[0] << std::endl;

			// build 3 strips, take best one each k
			size_t best_length = 0;
			size_t best_idx = 0;

			// init to gen 3 strips from face f
			HEMesh::HalfedgeHandle _start_hh;
			std::vector<HEMesh::FaceHandle> _faces[3];
			std::vector<HEMesh::FaceHandle> used_faces;
			HEMesh::HalfedgeHandle hh;
			HEMesh::FaceHandle fh;

			for (size_t i = 0; i < 3; ++i) {
				// build_strip(h[i], experiments[i], faces[i]);
				std::list<unsigned int> strip; // local
				_start_hh = h[i];

				strip.push_back(mesh.from_vertex_handle(_start_hh).idx());
				strip.push_back(mesh.to_vertex_handle(_start_hh).idx());

				// walk along the strip: 1st direction
				hh =
					mesh.prev_halfedge_handle(mesh.opposite_halfedge_handle(_start_hh));
				while (1) {
					// go right
					hh = mesh.next_halfedge_handle(hh);
					hh = mesh.opposite_halfedge_handle(hh);
					hh = mesh.next_halfedge_handle(hh);
					if (mesh.is_boundary(hh))
						break;
					fh = mesh.face_handle(hh);
					if (mesh.property(perFaceStripIdProperty, fh) != -1 ||
						mesh.property(faceUsedProperty, fh) == true)
						break;
					_faces[i].push_back(fh);
					// used(fh) = true;
					mesh.property(faceUsedProperty, fh) = true;
					strip.push_back(mesh.to_vertex_handle(hh).idx());

					// go left
					hh = mesh.opposite_halfedge_handle(hh);
					hh = mesh.next_halfedge_handle(hh);
					if (mesh.is_boundary(hh))
						break;
					fh = mesh.face_handle(hh);
					if (mesh.property(perFaceStripIdProperty, fh) != -1 ||
						mesh.property(faceUsedProperty, fh) == true)
						break;
					_faces[i].push_back(fh);
					mesh.property(faceUsedProperty, fh) = true;
					strip.push_back(mesh.to_vertex_handle(hh).idx());
				}

				// walk along the strip: 2nd direction
				bool flip(false);
				hh = mesh.prev_halfedge_handle(_start_hh);
				while (1) {
					// go right
					hh = mesh.next_halfedge_handle(hh);
					hh = mesh.opposite_halfedge_handle(hh);
					hh = mesh.next_halfedge_handle(hh);
					if (mesh.is_boundary(hh))
						break;
					fh = mesh.face_handle(hh);
					if (mesh.property(perFaceStripIdProperty, fh) != -1 ||
						mesh.property(faceUsedProperty, fh) == true)
						break;
					_faces[i].push_back(fh);
					mesh.property(faceUsedProperty, fh) = true;
					strip.push_front(mesh.to_vertex_handle(hh).idx());
					flip = true;

					// go left
					hh = mesh.opposite_halfedge_handle(hh);
					hh = mesh.next_halfedge_handle(hh);
					if (mesh.is_boundary(hh))
						break;
					fh = mesh.face_handle(hh);
					if (mesh.property(perFaceStripIdProperty, fh) != -1 ||
						mesh.property(faceUsedProperty, fh) == true)
						break;
					_faces[i].push_back(fh);
					mesh.property(faceUsedProperty, fh) = true;
					strip.push_front(mesh.to_vertex_handle(hh).idx());
					flip = false;
				}
				if (flip)
					strip.push_front(strip.front());

				const size_t length = strip.size();
				if (length > best_length) {
					best_length = length;
					best_idx = i;
				}

				for (auto f : mesh.faces()) {
					mesh.property(faceUsedProperty, f) = false;
				}
			}

			std::cout << "best index found for one face f random :" << best_idx
				<< std::endl;
			nStrips++;
			for (int i = 0; i < _faces[best_idx].size(); i++) {
				mesh.property(perFaceStripIdProperty, _faces[best_idx][i]) = nStrips;
				sampleset.remove(_faces[best_idx][i]);
				std::cout << "the length of a strip" << _faces[best_idx].size()
					<< std::endl;
			}
		}
	}

	return nStrips;
}

// operation 
void zr_add_noise(HEMesh& m) {
	std::mt19937 rnd;
	std::normal_distribution<float> dist;

	for (auto v : m.vertices()) {
		OpenMesh::Vec3f n;
		m.calc_vertex_normal_correct(v, n); // normal scales with area
		float areaScale = n.norm();
		float lengthScale = sqrt(areaScale);
		n = lengthScale / areaScale * n;

		m.point(v) += 0.1f * dist(rnd) * n;
	}
}
void zr_smooth_uniform_laplacian(HEMesh& m, float lambda, unsigned int iterations) {
	/*Task 2.2.4*/
	// this vector stores the computed centers of gravity
	std::vector<HEMesh::Point> centers_of_p;
	std::vector<HEMesh::Point>::iterator cog_it;
	centers_of_p.reserve(m.n_vertices());

	HEMesh::VertexIter v_it, v_end(m.vertices_end());
	HEMesh::VertexVertexIter vv_it;
	HEMesh::Point center_of_p;
	HEMesh::Scalar valence;

	OpenMesh::VPropHandleT<float> weights;
	m.add_property(weights);
	for (v_it = m.vertices_begin(); v_it != v_end; ++v_it) {
		m.property(weights, *v_it) = 0;
	}

	unsigned int i;
	for (i = 0; i < iterations; ++i) {
		centers_of_p.clear();
		/**
			  simple approch: set the weights to 1 (should from 0 to 1)
		*/
		for (v_it = m.vertices_begin(); v_it != v_end; ++v_it) {
			m.property(weights, *v_it) = 1;
		}

		// calculate centers with weights
		// this will not changed
		for (v_it = m.vertices_begin(); v_it != v_end; ++v_it) {
			center_of_p[0] = center_of_p[1] = center_of_p[2] = 0; // c[x] c[y] c[z]=0
			float deno = 0;
			for (vv_it = m.vv_iter(*v_it); vv_it.is_valid(); ++vv_it) {
				center_of_p += m.property(weights, *vv_it) * m.point(*vv_it);
				deno += m.property(weights, *vv_it);
				// std::cout << "deno:" << deno << std::endl;
			}
			if (deno == 0)
				deno = 1;
			centers_of_p.push_back(center_of_p / deno);
		}

		for (v_it = m.vertices_begin(), cog_it = centers_of_p.begin();
			v_it != v_end; ++v_it, ++cog_it) {
			*cog_it = *cog_it * lambda + m.point(*v_it) * (1 - lambda);
			m.set_point(*v_it, *cog_it);
		}
	}

	// do some test here

	// m.request_face_status();
	// m.request_edge_status();
	// m.request_vertex_status();

	// HEMesh::EdgeHandle eh = *m.edges_begin();
	// m.delete_edge(eh, true);
	// m.garbage_collection();

	int num = 0;
	for (HEMesh::EdgeIter e_it = m.edges_begin(); e_it != m.edges_end(); ++e_it) {
		++num;
	}
	std::cout << "edge num:" << num << std::endl;
}
void zr_smooth_cot_laplacian(HEMesh& m, float lambda, unsigned int iterations) {
	/*Task 2.2.4*/
	// this vector stores the computed centers of gravity
	std::vector<HEMesh::Point> centers_of_p;
	std::vector<HEMesh::Point>::iterator cog_it;
	centers_of_p.reserve(m.n_vertices());

	HEMesh::VertexIter v_it, v_end(m.vertices_end());
	HEMesh::VertexVertexIter vv_it;
	HEMesh::Point center_of_p;
	HEMesh::Scalar valence;

	OpenMesh::VPropHandleT<float> weights;
	m.add_property(weights);
	for (v_it = m.vertices_begin(); v_it != v_end; ++v_it) {
		m.property(weights, *v_it) = 0;
	}

	unsigned int i;
	for (i = 0; i < iterations; ++i) {
		centers_of_p.clear();
		/**
			  simple approch: set the weights to 1 (should from 0 to 1)
		*/

		// for (v_it = m.vertices_begin(); v_it != v_end; ++v_it) {
		//   m.property(weights, *v_it) = 1;
		// }

		/**
			  approch: set weights based on valence ?
		*/
		// for (v_it = m.vertices_begin(); v_it != v_end; ++v_it) {
		//   for (vv_it = m.vv_iter(*v_it); vv_it.is_valid(); ++vv_it) {
		//     m.property(weights, *v_it) += 1.0;
		//   }
		// }
		// for (v_it = m.vertices_begin(); v_it != v_end; ++v_it) {
		//   m.property(weights, *v_it) = 1.0 / m.property(weights, *v_it);
		// }

		/**
			  approch: set weights based on cot
		*/
		HEMesh::HalfedgeHandle heh0, heh1, heh2;
		HEMesh::VertexHandle v0, v1;
		HEMesh::Point* p0, * p1, * p2;
		HEMesh::Normal d0, d1;
		HEMesh::Scalar weight, lb(-1.0), ub(1.0);
		// each time , two points
		for (HEMesh::EdgeIter e_it = m.edges_begin(); e_it != m.edges_end();
			++e_it) {
			weight = 0.0;

			heh0 = m.halfedge_handle(*e_it, 0);
			v0 = m.to_vertex_handle(heh0);
			p0 = &m.point(v0);

			heh1 = m.halfedge_handle(*e_it, 1);
			v1 = m.to_vertex_handle(heh1);
			p1 = &m.point(v1);

			heh2 = m.next_halfedge_handle(heh0);
			p2 = &m.point(m.to_vertex_handle(heh2));
			d0 = (*p0 - *p2);
			d0.normalize();
			d1 = (*p1 - *p2);
			d1.normalize();
			weight += 1.0 / tan(acos(std::max(lb, std::min(ub, (d0 | d1)))));

			heh2 = m.next_halfedge_handle(heh1);
			p2 = &m.point(m.to_vertex_handle(heh2));
			d0 = (*p0 - *p2);
			d0.normalize();
			d1 = (*p1 - *p2);
			d1.normalize();
			weight += 1.0 / tan(acos(std::max(lb, std::min(ub, (d0 | d1)))));

			weight = weight / 2.0;
			m.property(weights, v0) += weight;
			m.property(weights, v1) += weight;
		}
		for (v_it = m.vertices_begin(); v_it != v_end; ++v_it) {
			if (m.property(weights, *v_it))
				m.property(weights, *v_it) = 1.0 / m.property(weights, *v_it);
			// std::cout << "weights each point:" << m.property(weights, *v_it)
			//           << std::endl;
		}

		// calculate centers with weights
		// this will not changed
		for (v_it = m.vertices_begin(); v_it != v_end; ++v_it) {
			center_of_p[0] = center_of_p[1] = center_of_p[2] = 0; // c[x] c[y] c[z]=0
			float deno = 0;
			for (vv_it = m.vv_iter(*v_it); vv_it.is_valid(); ++vv_it) {
				center_of_p += m.property(weights, *vv_it) * m.point(*vv_it);
				deno += m.property(weights, *vv_it);
				// std::cout << "deno:" << deno << std::endl;
			}
			if (deno == 0)
				deno = 1;
			centers_of_p.push_back(center_of_p / deno);
		}

		for (v_it = m.vertices_begin(), cog_it = centers_of_p.begin();
			v_it != v_end; ++v_it, ++cog_it) {
			*cog_it = *cog_it * lambda + m.point(*v_it) * (1 - lambda);
			m.set_point(*v_it, *cog_it);
		}
	}
}

// checking 
void zr_dualcalcu(HEMesh& m) {
	// Task

	HEMesh dual;

	OpenMesh::FPropHandleT<HEMesh::VertexHandle> centers;
	m.add_property(centers);

	// loop all faces, add vertex to dual
	for (HEMesh::ConstFaceIter fit = m.faces_begin(); fit != m.faces_end();
		++fit) {
		HEMesh::Point centerPoint(0, 0, 0);
		HEMesh::Scalar degree = 0.0;
		for (HEMesh::ConstFaceVertexIter vit = m.cfv_iter(*fit); vit.is_valid();
			++vit, ++degree)
			centerPoint += m.point(*vit);
		assert(degree != 0);
		centerPoint /= degree;
		m.property(centers, *fit) = dual.add_vertex(centerPoint);
	}

	// loop all vertex in m, add a face in the dual
	std::vector<HEMesh::VertexHandle> face_vhandles;
	for (HEMesh::ConstVertexIter vit = m.vertices_begin();
		vit != m.vertices_end(); ++vit) {
		if (!m.is_boundary(*vit)) {
			face_vhandles.clear();
			for (HEMesh::ConstVertexFaceIter fit = m.cvf_iter(*vit); fit.is_valid();
				++fit)
				face_vhandles.push_back(m.property(centers, *fit));
			dual.add_face(face_vhandles);
		}
	}

	m.remove_property(centers);

	m.clear();
	m = dual;
}
bool zr_is_dege_flip_ok(HEMesh& m, HEMesh::EdgeHandle eh) {

	if (m.is_boundary(eh))
		return false;
	else {
		HEMesh::HalfedgeHandle heh0, heh1, heh2, heh3;
		HEMesh::VertexHandle v0, v1, v2, v3;
		HEMesh::Point* p0, * p1, * p2;
		heh0 = m.halfedge_handle(eh, 0);
		v0 = m.to_vertex_handle(heh0);
		p0 = &m.point(v0);

		heh1 = m.halfedge_handle(eh, 1);
		v1 = m.to_vertex_handle(heh1);
		p1 = &m.point(v1);

		heh2 = m.next_halfedge_handle(heh0);
		v2 = m.to_vertex_handle(heh2);

		heh3 = m.next_halfedge_handle(heh1);
		v3 = m.to_vertex_handle(heh3);

		for (HEMesh::VertexEdgeIter eitv2 = m.ve_iter(v2); eitv2.is_valid();
			++eitv2) {
			for (HEMesh::VertexEdgeIter eitv3 = m.ve_iter(v3); eitv3.is_valid();
				++eitv3) {
				if (*eitv2 == *eitv3)
					return false;
			}
		}
	}

	return true;
}
void zr_flipcheck(HEMesh& m) {
	// Task
	m.request_face_status();
	m.request_edge_status();
	m.request_vertex_status();

	for (HEMesh::EdgeIter e_it = m.edges_begin(); e_it != m.edges_end(); ++e_it) {
		if (zr_is_dege_flip_ok(m, *e_it)) {
			m.delete_edge(*e_it, false);
		}
	}

	m.garbage_collection();
}

// computing 
float zr_calcu_surfacearea(const HEMesh& m) {
	float area = 0;
	/* Task 2.2.2 */
	int face = 0;
	int points_repeat = 0;
	for (HEMesh::VertexIter v_it = m.vertices_begin(); v_it != m.vertices_end();
		++v_it) {
		for (HEMesh::FaceHandle vf_it : m.vf_range(v_it)) { // get into a face
			float polygonarea = 0;
			face += 1;

			std::vector<HEMesh::Point> points;
			for (HEMesh::VertexHandle p :
			m.fv_range(vf_it)) { // add the points in a surface into a list
				points.push_back(m.point(p));
			}
			points_repeat = points.size();
			std::cout << "points in current face :" << points_repeat << std::endl;
			// when it is triangle mesh , the size should be 3 for all faces
			for (int j = 0; j < points_repeat - 2; j++) {
				polygonarea +=
					((points[j + 1] - points[0]) % (points[j + 2] - points[0])).norm() /
					2;
			}
			area += polygonarea / points_repeat; // we calcu area points_repeat times
		}
	}
	return area;
}
float zr_calcu_volume(const HEMesh& m) {
	float vol = 0;
	/*Task 2.2.2*/
	// std::cout << "Volume calculation is not implemented." << std::endl;
	int face = 0;
	int points_repeat = 0;
	for (HEMesh::VertexIter v_it = m.vertices_begin(); v_it != m.vertices_end();
		++v_it) {
		for (HEMesh::FaceHandle vf_it : m.vf_range(v_it)) { // get into a face
			float polygonvol = 0;
			face += 1;

			std::vector<HEMesh::Point> points;
			for (HEMesh::VertexHandle p : m.fv_range(vf_it)) { // add the points in a surface into a list
				points.push_back(m.point(p));
			}
			points_repeat = points.size();
			std::cout << "points in current face :" << points_repeat << std::endl;
			for (int j = 0; j < points_repeat - 2; j++) {
				polygonvol += (points[0] | ((points[j + 1] - points[0]) %
					(points[j + 2] - points[0]))) /
					6;
			}
			vol += polygonvol / points_repeat; // we calcu area points_repeat times
		}
	}
	std::cout << "faces in model:" << face / points_repeat << std::endl;

	return vol;
}
bool zr_compute_ray_tri_intersection(OpenMesh::Vec3f ori, OpenMesh::Vec3f raydir, OpenMesh::Vec3f plane_normal) {
	float t, s;
	t = dot(plane_normal, ori);
	s = dot(plane_normal, raydir);
	if(s == 0)
		return false;
	t = 1 - t;
	t /= s;
	return (t >= 0);
}
void zr_compute_mesh_ray_intersected_tri(HEMesh& m, OpenMesh::FPropHandleT<int> perFaceShellIndex, OpenMesh::Vec3f ori, OpenMesh::Vec3f intersectionP) {
	for (auto f : m.faces()) {
		m.property(perFaceShellIndex, f) = 1;
	};
	std::vector<HEMesh::Point> points;
	OpenMesh::Vec3f n;
	for (auto f : m.faces()) {
		points.clear();
		n = m.calc_face_normal(f); 
		n.normalize();
		/*if (zr_compute_ray_tri_intersection(ori, raydir, n)) {
			m.property(perFaceShellIndex, f) = 2;
			break;
		}*/
		for (auto p : m.fv_range(f)) {
			points.push_back(m.point(p));
		}

		//this could not work other than tri. meshes, but can be improved
		OpenMesh::Vec3f A;
		OpenMesh::Vec3f B;
		OpenMesh::Vec3f C;
		OpenMesh::Vec3f Q = intersectionP;
		//intersectionP

		for (int i = 2; i < points.size(); i++) {
			A = points.at(i-2);
			B = points.at(i-1);
			C = points.at(i);
			if ((dot(cross(B - A, Q - A), n) >= 0) && (dot(cross(C - B, Q - B), n) >= 0) && (dot(cross(A - C, Q - C), n) >= 0)) {
				// inside!
				m.property(perFaceShellIndex, f) = 2;
				//break;
			}
		}

	};
}

#endif
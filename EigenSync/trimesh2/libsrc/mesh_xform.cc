#include "TriMesh.h"
#include "TriMesh_algo.h"
#include "XForm.h"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cctype>
using namespace std;
using namespace trimesh;

#define ATOF(x) ((float) atof(x))


void apply_xform(TriMesh *mesh, const char *xffilename)
{
	xform xf;
	if (!xf.read(xffilename))
		fprintf(stderr, "Couldn't open %s\n", xffilename);
	else
		apply_xform(mesh, xf);
}

void apply_xform2(TriMesh *mesh, const xform &xf)
{
	int nv = mesh->vertices.size();

#pragma omp parallel for
	for (int i = 0; i < nv; i++)
		mesh->vertices[i] = xf * mesh->vertices[i];

	if (!mesh->normals.empty()) {
		xform nxf = norm_xf(xf);
#pragma omp parallel for
		for (int i = 0; i < nv; i++) {
			mesh->normals[i] = nxf * mesh->normals[i];
			normalize(mesh->normals[i]);
		}
	}

	if (mesh->bbox.valid) {
		mesh->bbox.valid = false;
		mesh->need_bbox();
	}
	if (mesh->bsphere.valid) {
		mesh->bsphere.valid = false;
		mesh->need_bsphere();
	}
}

int main(int argc, char const *argv[])
{
	
	return 0;
}
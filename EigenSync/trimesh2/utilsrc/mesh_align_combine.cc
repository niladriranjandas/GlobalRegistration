/*
Szymon Rusinkiewicz
Princeton University

mesh_align.cc
Minimal interface to ICP: register two meshes given an initial guess
for their alignment.
*/
#include <fstream>
#include <iostream>
#include "TriMesh.h"
#include "TriMesh_algo.h"
#include "ICP.h"
#include "ICP.cc"
#include <cstdio>
#include <cstdlib>
#include <unistd.h>
#include <vector>
#include <string>
using namespace std;
using namespace trimesh;
 //vector<PtPair> pairs;






void usage(const char *myname)
{
	fprintf(stderr, "Usage: %s [-options] mesh1.ply mesh2.ply\n", myname);
	fprintf(stderr, "Reads transforms in mesh1.xf and mesh2.xf, updates the latter\n");
	fprintf(stderr, "Options:\n");
	fprintf(stderr, "	-a		Align using affine xform\n");
	fprintf(stderr, "	-r		Align using rigid-body transform (default)\n");
	fprintf(stderr, "	-s		Align using rigid + isotropic scale\n");
	fprintf(stderr, "	-v		Verbose\n");
	fprintf(stderr, "	-b		Bulk mode: overlap checking, write to mesh1--mesh2.xf\n");
	exit(1);
}

int main(int argc, char *argv[])
{
	int verbose = 0;
	bool do_scale = false;
	bool do_affine = false;
	bool bulkmode = false;

	int c;
	while ((c = getopt(argc, argv, "harsvb")) != EOF) {
		switch (c) {
			case 'a': do_affine = true; do_scale = false; break;
			case 'r': do_affine = do_scale = false; break;
			case 's': do_scale = true; do_affine = false; break;
			case 'v': verbose = 2; break;
			case 'b': bulkmode = true; break;
			default: usage(argv[0]);
		}
	}

	TriMesh::set_verbose(verbose);

	if (argc - optind < 2)
		usage(argv[0]);
	const char *filename1 = argv[optind], *filename2 = argv[optind+1];

	TriMesh *mesh1 = TriMesh::read(filename1);
	if (!mesh1)
		usage(argv[0]);
	TriMesh *mesh2 = TriMesh::read(filename2);
	if (!mesh2)
		usage(argv[0]);

	xform xf1;
	string xffilename1 = xfname(filename1);
	xf1.read(xffilename1);

	xform xf2;
	string xffilename2 = xfname(filename2);
	xf2.read(xffilename2);

	KDtree *kd1 = new KDtree(mesh1->vertices);
	KDtree *kd2 = new KDtree(mesh2->vertices);
	vector<float> weights1, weights2;

	if (bulkmode) {
		float area1 = mesh1->stat(TriMesh::STAT_TOTAL, TriMesh::STAT_FACEAREA);
		float area2 = mesh2->stat(TriMesh::STAT_TOTAL, TriMesh::STAT_FACEAREA);
		float overlap_area, overlap_dist;
		find_overlap(mesh1, mesh2, xf1, xf2, kd1, kd2,
			overlap_area, overlap_dist);
		float frac_overlap = overlap_area / min(area1, area2);
		if (frac_overlap < 0.1f) {
			TriMesh::eprintf("Insufficient overlap\n");
			exit(1);
		} else {
			TriMesh::dprintf("%.1f%% overlap\n",
				frac_overlap * 100.0);
		}
	}

	float err = ICP(mesh1, mesh2, xf1, xf2, kd1, kd2, weights1, weights2,
		verbose, do_scale, do_affine);
	if (err >= 0.0f){
		err = ICP(mesh1, mesh2, xf1, xf2, kd1, kd2, weights1, weights2,
			verbose, do_scale, do_affine);
	
	}

	if (err < 0.0f) {
		TriMesh::eprintf("0\n");
		exit(1);
	}

	TriMesh::eprintf("%f\n", err);




		size_t nv1 = mesh1->vertices.size(), nv2 = mesh2->vertices.size();

	vector<float> sampcdf1(nv1), sampcdf2(nv2);
	for (size_t i = 0; i < nv1-1; i++)
		sampcdf1[i] = (float) (i+1) / nv1;
	sampcdf1[nv1-1] = 1.0f;
	for (size_t i = 0; i < nv2-1; i++)
		sampcdf2[i] = (float) (i+1) / nv2;
	sampcdf2[nv2-1] = 1.0f;
    
 vector<PtPair> pairs;
    float incr =4.0/500;
   
    select_and_match(mesh1, mesh2, xf1, xf2, kd2, sampcdf1, incr,
			 0.0, verbose, pairs, false);
	select_and_match(mesh2, mesh1, xf2, xf1, kd1, sampcdf2, incr,
			 0.0, verbose, pairs, true);

    size_t np = pairs.size();
 
		//TriMesh::eprintf("Generated %lu pairs .\n",
		//	(unsigned long) np);
	

	// Reject pairs with distance > 3 sigma
	float thresh = 19.772984f * median_dist2(pairs);

///		TriMesh::eprintf("Rejecting pairs > %f\n", sqrt(thresh));
	size_t next = 0;
	for (size_t i = 0; i < np; i++) {
		if (dist2(pairs[i].p1, pairs[i].p2) <= thresh)
			pairs[next++] = pairs[i];
	}
	pairs.erase(pairs.begin() + next, pairs.end());
	np = pairs.size();
	// for (size_t i = 0; i < nv1; i++)
	// cout<<mesh1->vertices[i]<<endl;

   




	if (bulkmode) {
		string xffilename12 = filename1;
		size_t dot = xffilename12.rfind(".", xffilename12.length());
		if (dot != string::npos)
			xffilename12.erase(dot);
		xffilename12 += string("_") + replace_ext(filename2, "xf");
		xform xf12 = inv(xf2) * xf1;
		xf12.write(xffilename12);
	} else {
		string xffilename12 = filename1;
		size_t dot = xffilename12.rfind(".", xffilename12.length());
		if (dot != string::npos)
			xffilename12.erase(dot);
		xffilename12 += string("_") + replace_ext(filename2, "xf");
		xf2.write(xffilename12);


        string xffilename123 = filename1;
        size_t dot2 = xffilename123.rfind(".", xffilename123.length());
        if (dot != string::npos)
			xffilename123.erase(dot2);		xffilename123 += string("_");

        xffilename123 += string("_") + replace_ext(filename2, "pair");
		    ofstream myfile;
            myfile.open (xffilename123.c_str());
            	np = pairs.size();
				for (size_t i = 0; i < np; i++)
          			 myfile<<pairs[i].p1<<", "<<pairs[i].p2<<endl;








	}
	exit(0);
}


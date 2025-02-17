#include <iostream>
#include <fstream>
//#include <iostream>
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






int main(int argc, char *argv[])
{
	int verbose = 0;
	bool do_scale = false;
	bool do_affine = false;
	bool bulkmode = false;
    //ofstream myfile;
	int c;
	TriMesh::set_verbose(verbose);
    const char *filename1 = argv[1];
	TriMesh *mesh1 = TriMesh::read(filename1);
	
	xform xf1;
	string xffilename1 = xfname(filename1);
	xf1.read(xffilename1);

	size_t nv1 = mesh1->vertices.size();
    string xffilename12 = filename1;
	size_t dot = xffilename12.rfind(".", xffilename12.length());
	if (dot != string::npos)
			xffilename12.erase(dot);
	

   xffilename12 =replace_ext(filename1, "vtx");

//  cout<<xffilename12<<endl;

  std::ofstream myfile;//(xffilename12);
  myfile.open("hello.vtx");
  	for (size_t i = 0; i < nv1; i++)
	myfile<<mesh1->vertices[i]<<endl;
 myfile.close();

   return 0;
}
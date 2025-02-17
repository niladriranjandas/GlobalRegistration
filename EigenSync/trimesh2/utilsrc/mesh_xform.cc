#include "TriMesh.h"
#include "TriMesh_algo.h"
#include "XForm.h"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cctype>
#include <ctime>           

#include <boost/random/linear_congruential.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/generator_iterator.hpp>
using namespace std;
using namespace trimesh;

typedef boost::minstd_rand base_generator_type;
base_generator_type generator(912);
float uni=0;
int itm=0;




// Transform the mesh by the given matrix
void apply_xform1(TriMesh *mesh, const xform &xf)
{cout<<"hello3";
     
     boost::uniform_real<> uni_dist(-1*uni,uni);
     boost::variate_generator<base_generator_type&, boost::uniform_real<> > uni(generator, uni_dist);  
	
     int nv = mesh->vertices.size();
     float noise =0;
#pragma omp parallel for
	for (int i = 0; i < nv; i++)
		{mesh->vertices[i] = xf * mesh->vertices[i];
                 if(itm==0)
                    noise=0;
                 else  
                    noise=uni();
                 mesh->vertices[i][0]=mesh->vertices[i][0]+noise;
                 mesh->vertices[i][1]=mesh->vertices[i][1]+noise;
                 mesh->vertices[i][2]=mesh->vertices[i][2]+noise;   

                    }

	if (!mesh->normals.empty()) {
		xform nxf = norm_xf(xf);
#pragma omp parallel for
		for (int i = 0; i < nv; i++) {
			mesh->normals[i] = nxf * mesh->normals[i];
			normalize(mesh->normals[i]);
                 if(itm==0)
                    noise=0;
                 else  
                    noise=uni();
                          
                           
                             
                       mesh->normals[i][0]=mesh->normals[i][0]+noise;
                       mesh->normals[i][1]=mesh->normals[i][1]+noise;
                       mesh->normals[i][2]=mesh->normals[i][2]+noise;
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





void apply_xform(TriMesh *mesh, const char *xffilename)
{cout<<"hello2";
	xform xf;
	if (!xf.read(xffilename))
		fprintf(stderr, "Couldn't open %s\n", xffilename);
	else
		apply_xform1(mesh, xf);
}



int main(int argc, char const *argv[])
{

	
 	const char  *filename1 = argv[1];
 	const char  *filename2 = argv[2];
    const char *outfile = argv[3];
    uni= atof(argv[4]);
    itm=atoi(argv[5]);
   		 if(itm==0)
                    uni=1;
                     

    cout<<"hello"<<uni<<endl;
    TriMesh *mesh1 = TriMesh::read(filename1);
	apply_xform(mesh1,filename2);
    mesh1->write(outfile);
	return 0;
}

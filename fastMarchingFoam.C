
/*
Notes:
- required to comment line 1164 ("assert(p->start() < p->stop());") out in geodesic/geodesic_algorithm_exact.h
  this is also done in the AIMS 4.4 library (http://brainvisa.info/downloadpage.html)
*/

#include <vector>
#include <string>
#include <cstring>
#include <algorithm>
#include <memory>
#include "geodesic/geodesic_algorithm_exact.h"
#include "geodesic/geodesic_algorithm_dijkstra.h"
#include "fvCFD.H"
#include "calculatedFvPatchField.H"


#define FM_HUGE 1e20


class GeodesicFoamMesh
{
protected:
	// maps global face indices to (geodesic) face center vertex indices
	HashTable<label, label> face_to_Cf_vertex_map;

	// mesh and algorithm objects
	const Foam::fvMesh* pmesh;
	geodesic::Mesh* pgeo_mesh;
	geodesic::GeodesicAlgorithmBase* pgeo_algorithm;

	// vertices and triangles
	std::vector<vector> vertices;
	std::vector<std::size_t> triangles;

public:
	GeodesicFoamMesh(const Foam::fvMesh& mesh)
	{
		pmesh = &mesh;

		// get the boundary mesh
		const polyBoundaryMesh& bmesh = mesh.boundaryMesh();

		// get mesh points (vertices) and faces
		const pointField& points = mesh.points();
		const faceList& faces = mesh.faces();

		// hash table for testing if a point was already added to vertices vector
		HashTable<label, label> point_to_vertex_map;

		// maps the global face index to the index of the face center vertex in the vertices vector
		HashTable<label, label> face_to_vertex_map;

		//std::cout << Pstream::myProcNo() << " points " << points.size() << std::endl;

		// iterate over every patch of the boundary mesh
		forAll(bmesh, patchI)
		{
			// get the patch
			const polyPatch& patch = bmesh[patchI];

			// get the face centers
			const vectorField& Cf_patch = mesh.Cf().boundaryField()[patchI];

			// get the face normals
			const vectorField& Sf_patch = mesh.Sf().boundaryField()[patchI];

			//std::cout << Pstream::myProcNo() << " patch " << patch.name() << ": " << patch.size() << std::endl;

			// iterate over all faces of the patch
			forAll(patch, faceI)
			{
				// calculate the global face index (on this processor)
				label globalFaceI = patch.start() + faceI;

				// get the face
				const face& f = faces[globalFaceI];

				// get face center
				vector Cf = Cf_patch[faceI];

				vector Sf = Sf_patch[faceI];
				vector n = Sf / mag(Sf);

				//std::cout << Pstream::myProcNo() << " face " << faceI << " " << globalFaceI << std::endl;

				// get the first point
				vector p0 = points[f[0]];

				// map for angle of face point (relative to center and first point) to global face point index
				HashTable<label, scalar> angle_to_point_map;

				// iterate over face points
				forAll(f, pointI)
				{
					// get the global point index (on this processor)
					label globalPointI = f[pointI];

					// get the point
					vector p = points[globalPointI];

					// compute the angle between (p0 - Cf) and (p - Cf)
					vector Cfp0 = p0 - Cf; Cfp0 /= mag(Cfp0);
					vector Cfp = p - Cf; Cfp /= mag(Cfp);
					scalar cosAngle = Cfp0 & Cfp;
					scalar sinAngle = (Cfp0 ^ Cfp) & n;
					scalar angle = std::atan2(sinAngle, cosAngle);
					if (angle < 0) angle += 2*M_PI;
					angle_to_point_map.insert(angle, globalPointI);

					// add point to liCfp0st of vertices
					if (!point_to_vertex_map.found(globalPointI)) {
						// make sure to add this point not again
						point_to_vertex_map.insert(globalPointI, vertices.size());
						// add point to vertex list
						vertices.push_back(points[globalPointI]);
					}
				}

				// triangulate the face
				// we sort the face points by angle, since we do not know if they are sorted correctly
				List<scalar> angles = angle_to_point_map.sortedToc();

				// add face center to vertex list
				label Cf_vertex_index = vertices.size();
				vertices.push_back(Cf);
				face_to_Cf_vertex_map.insert(globalFaceI, Cf_vertex_index);

				forAll(f, angleI)
				{
					// get the triangle vertices
					label v0 = Cf_vertex_index;
					label v1 = point_to_vertex_map[angle_to_point_map[angles[angleI]]];
					label v2 = point_to_vertex_map[angle_to_point_map[angles[(angleI+1) % angles.size()]]];

					//std::cout << Pstream::myProcNo() << " triangle " << v0 << " " << v1 << " " << v2 << std::endl;

					// add the triangle
					triangles.push_back(v0);
					triangles.push_back(v1);
					triangles.push_back(v2);
				}

			}
		}

		// convert vetrices and triangles to geodesic format
		std::vector<double> geo_vertices(vertices.size()*3);
		std::vector<unsigned> geo_triangles(triangles.size());

		for (std::size_t i = 0; i < vertices.size(); i++) {
			geo_vertices[3*i + 0] = (double) vertices[i].x();
			geo_vertices[3*i + 1] = (double) vertices[i].y();
			geo_vertices[3*i + 2] = (double) vertices[i].z();
		}

		for (std::size_t i = 0; i < triangles.size(); i++) {
			geo_triangles[i] = (unsigned) triangles[i];
		}

		// create internal mesh data structure including edges
		pgeo_mesh = new geodesic::Mesh();
		pgeo_mesh->initialize_mesh_data(geo_vertices, geo_triangles);		

		// create exact algorithm for the mesh
		pgeo_algorithm = new geodesic::GeodesicAlgorithmExact(pgeo_mesh);
	}

	~GeodesicFoamMesh()
	{
		delete pgeo_mesh;
		delete pgeo_algorithm;
	}

	label closestVertex(const vector& v)
	{
		label imin = 0;
		double dmin = Foam::mag(vertices[0] - v);

		for (std::size_t i = 1; i < vertices.size(); i++) {
			vector dv = vertices[i] - v;
			scalar d = Foam::mag(dv);
			if (d < dmin) {
				dmin = d;
				imin = i;
			}
		}

		return imin;
	}

	void calculateDistance(volScalarField& d, label vertexIndex, double maxDist = FM_HUGE)
	{
		// get the boundary mesh
		const polyBoundaryMesh& bmesh = pmesh->boundaryMesh();

		// get mesh faces
		const faceList& faces = pmesh->faces();

		// create source 
		// this is the point with zero distance
		// in general, there could be multiple sources, but now we have only one
		unsigned source_vertex_index = vertexIndex;
		geodesic::SurfacePoint source(&pgeo_mesh->vertices()[source_vertex_index]);
		std::vector<geodesic::SurfacePoint> all_sources(1,source);					

		// calculated distances to all vertices (on surface)
		pgeo_algorithm->propagate(all_sources, maxDist);

		// set default value
		d.boundaryField() = maxDist;

		// save distance in field d
		forAll(bmesh, patchI)
		{
			const polyPatch& patch = bmesh[patchI];

			forAll(patch, faceI)
			{
				label globalFaceI = patch.start() + faceI;
				const face& f = faces[globalFaceI];

				// get geodesic vertex index
				label vertexI = face_to_Cf_vertex_map[globalFaceI];

				// define geodesic surface point
				geodesic::SurfacePoint p(&pgeo_mesh->vertices()[vertexI]);		

				// find closets source and distance to the source
				double distance;
				unsigned best_source = pgeo_algorithm->best_source(p, distance);

				// store distance in field d
				d.boundaryField()[patchI][faceI] = (scalar)std::min(distance, maxDist);
			}
		}
	}

#if 0
	void closeSurface(const std::vector<geodesic:SurfacePoint>& path)
	{
		
	}

	void cutPath(const std::vector<geodesic:SurfacePoint>& path)
	{

		// TODO: get common edge of path[i] and path[i+1] vertex

		// determine triangle left of path and right of path
		// by checking sign of ((path[i+1] - path[i]) ^ (v - path[i])) & n
		// where n is the triangle face normal and v is the vertex opposite to the triangle edge

		// replace points of left triangle edge by duplicates
		// except first and last vertex, if path is not closed
		// except last vertex, if path is closed

		
	}
#endif

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
	Foam::argList args(argc, argv);

	// check arguments
	if (!args.checkRootCase()) {
		Foam::FatalError.exit();
	}

	// create runtime object
	Foam::Time runTime(Foam::Time::controlDictName, args);

	// load the mesh
	Foam::fvMesh mesh(
		Foam::IOobject
		(
			Foam::fvMesh::defaultRegion,
			runTime.timeName(),
			runTime,
			Foam::IOobject::MUST_READ
		)
	);

	std::cout << "FM_HUGE = " << FM_HUGE << std::endl;

	Info<< "Reading fastMarchingDict\n" << endl;

	// read dictionary containing the points to be processed
	Foam::IOdictionary fastMarchingDict(
		Foam::IOobject(
			"fastMarchingDict",
			runTime.constant(),
			mesh,
			Foam::IOobject::MUST_READ,
			Foam::IOobject::NO_WRITE
		)
	);

	// iterate over the points and process each point
	forAll(fastMarchingDict.toc(), key)
	{
		word name = fastMarchingDict.toc()[key];
		dictionary& dict = fastMarchingDict.subDict(name);

		Info << "Processing fast marching point: " << name << endl;

		// read position
		vector center(dict.lookup("sourcePoint"));

		// read maximum propagation distance
		scalar maxDist(dict.lookupOrDefault<scalar>("maxDist", FM_HUGE));

		// calculated distance field (actually only the surface faces are used)
		volScalarField d
		(
			IOobject
			(
				name,
				runTime.timeName(),
				mesh,
				IOobject::NO_READ,
				IOobject::NO_WRITE
			),
			mesh, 0, "calculated"
		);

		// calculated distance
		GeodesicFoamMesh gfmesh(mesh);

		label vertex_index = gfmesh.closestVertex(center);
		gfmesh.calculateDistance(d, vertex_index, maxDist);

		// write field d
		d.write();
	}

	// program successful
	return 0;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //





// The following is all old stuff:


void getNeighboursRecursive(polyMesh& mesh, label cellI, int nmax, int levels, const labelHashSet& constraint, HashTable<scalar,label>& result);
void getNeighboursRecursive(polyMesh& mesh, label cellI, int nmax, int levels, const labelHashSet& constraint, HashTable<scalar,label>& result, label currentCellI, labelHashSet& visitedCells);




int xmain(int argc, char *argv[])
{
	Foam::argList args(argc, argv);

	// check arguments
	if (!args.checkRootCase()) {
		Foam::FatalError.exit();
	}

	// create runtime object
	Foam::Time runTime(Foam::Time::controlDictName, args);

	// load the mesh
	Foam::fvMesh mesh
	(
		Foam::IOobject
		(
			Foam::fvMesh::defaultRegion,
			runTime.timeName(),
			runTime,
			Foam::IOobject::MUST_READ
		)
	);


	volVectorField g
	(
		IOobject
		(
			"g",
			runTime.timeName(),
			mesh,
			IOobject::NO_READ,
			IOobject::NO_WRITE
		),
		mesh, vector(0,0,0), "zeroGradient"
	);

	volScalarField d
	(
		IOobject
		(
			"d",
			runTime.timeName(),
			mesh,
			IOobject::NO_READ,
			IOobject::NO_WRITE
		),
		mesh, 0, "calculated"
	);


/*
	d.boundaryField()[0][10] = 11;
	d.write();
	std::cout << "patchI " << d.boundaryField()[0].size() << " " << d.boundaryField()[0].type() << std::endl;
	return 0;
*/

	// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

	const polyBoundaryMesh& bmesh = mesh.boundaryMesh();
	const pointField& points = mesh.points();
	const faceList& faces = mesh.faces();
	HashTable<label, label> boundary_face_hash;
	HashTable<label, label> boundary_point_hash;
	labelList boundary_faces;
	labelList boundary_points;

	forAll(bmesh, patchI)
	{
		const polyPatch& patch = bmesh[patchI];

		// std::cout << "patch " << patch.name() << ": " << patch.size() << std::endl;

		forAll(patch, faceI)
		{
			label globalFaceI = patch.start() + faceI;
			// std::cout << "face " << faceI << " " << globalFaceI << std::endl;
			const face& f = faces[globalFaceI];

			if (!boundary_face_hash.found(globalFaceI)) {
				boundary_faces.append(globalFaceI);
				boundary_face_hash.insert(globalFaceI,  boundary_face_hash.size());
			}
		
			forAll(f, pointI)
			{
				label globalPointI = f[pointI];
				//std::cout << "point " << pointI << " " << globalPointI << std::endl;

				if (!boundary_point_hash.found(globalPointI)) {
					boundary_points.append(globalPointI);
					boundary_point_hash.insert(globalPointI,  boundary_point_hash.size());
				}
			}
		}
	}

	std::vector<double> vertices;
	std::vector<unsigned> triangles;

	// get vertices
	for (int i = 0; i < boundary_points.size(); i++)
	{
		const vector& p = points[boundary_points[i]];
		vertices.push_back(p.x());
		vertices.push_back(p.y());
		vertices.push_back(p.z());
	}

	// get faces
	for (int i = 0; i < boundary_faces.size(); i++)
	{
		const face& f = faces[boundary_faces[i]];
		
		int n = 0;
		forAll(f, pointI)
		{
			label globalPointI = f[pointI];

			if (n >= 3) {
				int m = triangles.size();
				triangles.push_back(triangles[m-n]);
				triangles.push_back(triangles[m-1]);
				triangles.push_back(boundary_point_hash[globalPointI]);
				n += 3;
			}
			else {
				triangles.push_back(boundary_point_hash[globalPointI]);
				n++;
			}
		}

		if (n < 3) {
			std::cout << "bad face " << i << std::endl;
			return -1;
		}
	}

	std::cout << "ntriangeles " << (triangles.size()/3) << std::endl;
	std::cout << "nvertices " << (vertices.size()/3) << std::endl;

#if 0
	for (size_t i = 0; i < vertices.size(); i++) {
		if ((i % 3) == 0) std::cout << i << ": ";
		std::cout << vertices[i] << " ";
		if (((i+1) % 3) == 0) std::cout << std::endl;
	}


	for (size_t i = 0; i < triangles.size(); i++) {
		if ((i % 3) == 0) std::cout << i << ": ";
		std::cout << triangles[i] << " ";
		if (((i+1) % 3) == 0) std::cout << std::endl;
	}
#endif

	geodesic::Mesh geo_mesh;
	geo_mesh.initialize_mesh_data(vertices, triangles);		//create internal mesh data structure including edges

	geodesic::GeodesicAlgorithmDijkstra algorithm(&geo_mesh);	//create exact algorithm for the mesh

	unsigned source_vertex_index = 0;

	geodesic::SurfacePoint source(&geo_mesh.vertices()[source_vertex_index]);		//create source 
	std::vector<geodesic::SurfacePoint> all_sources(1,source);					//in general, there could be multiple sources, but now we have only one

	if(false)	//target vertex specified, compute single path
	{
		unsigned target_vertex_index = atol(argv[3]);
		geodesic::SurfacePoint target(&geo_mesh.vertices()[target_vertex_index]);		//create source 

		std::vector<geodesic::SurfacePoint> path;	//geodesic path is a sequence of SurfacePoints

		bool const lazy_people_flag = false;		//there are two ways to do exactly the same
		if(lazy_people_flag)
		{
			algorithm.geodesic(source, target, path); //find a single source-target path
		}
		else		//doing the same thing explicitly for educational reasons
		{
			double const distance_limit = FM_HUGE;			// no limit for propagation
			std::vector<geodesic::SurfacePoint> stop_points(1, target);	//stop propagation when the target is covered
			algorithm.propagate(all_sources, distance_limit, &stop_points);	//"propagate(all_sources)" is also fine, but take more time because covers the whole mesh

			algorithm.trace_back(target, path);		//trace back a single path 
		}
		
		print_info_about_path(path);
		for(unsigned i = 0; i<path.size(); ++i)
		{
			geodesic::SurfacePoint& s = path[i];
			std::cout << s.x() << "\t" << s.y() << "\t" << s.z() << std::endl;
		}
	}
	else		//target vertex is not specified, print distances to all vertices
	{
		algorithm.propagate(all_sources);	//cover the whole mesh

#if 0
		for(unsigned i=0; i<geo_mesh.vertices().size(); ++i)
		{
			geodesic::SurfacePoint p(&geo_mesh.vertices()[i]);		

			double distance;
			unsigned best_source = algorithm.best_source(p,distance);		//for a given surface point, find closets source and distance to this source

			std::cout << distance << " ";		//print geodesic distance for every vertex
		}
		std::cout << std::endl;
#endif

		forAll(bmesh, patchI)
		{

			const polyPatch& patch = bmesh[patchI];
			const vectorField& cf = mesh.Cf().boundaryField()[patchI];


			forAll(patch, faceI)
			{
				label globalFaceI = patch.start() + faceI;
				// std::cout << "face " << faceI << " " << globalFaceI << std::endl;
				const face& f = faces[globalFaceI];
				
				d.boundaryField()[patchI][faceI] = 0.0;

//				Foam::fvPatchField<scalar>& bf = d.boundaryField()[patchI];
//				bf[faceI] = 0.0;

				vector c = cf[faceI]; 

//				std::cout << "face " << faceI << " " << c.x() << std::endl;

//				geodesic::SurfacePoint p(&geo_mesh.vertices()[i]);		
			}

		}
	}


	return 0;



	// this is the old code (not very accurate, but shows how to work with cells)

#if 1
	const vectorField& c =	mesh.cellCentres();

	labelHashSet frozenCells;
	labelHashSet newFront;

	frozenCells.insert(23415);

	// set inital distance value
	for (labelHashSet::iterator i = frozenCells.begin(); i != frozenCells.end(); i++) {
		d[i.key()] = 0.0;
		g[i.key()] = vector(0.0, 0.0, 0.0);
	}

	for(;;)
	{
		// recompute new front
		newFront.clear();
		for (labelHashSet::iterator i = frozenCells.begin(); i != frozenCells.end(); i++)
		{
			// get cell neighbours
			label cellI = i.key();
			const labelList& neighbours = mesh.cellCells()[cellI];

			// iterate over cell neighbours
			for (int j = 0; j < neighbours.size(); j++)
			{
				label neighbourCellI = neighbours[j];

				// skip already computed cells
				if (d[neighbourCellI] >= 0.0) {
					continue;
				}
				newFront.set(neighbourCellI);
			}
		}

		// check if we are done
		if (newFront.size() == 0) break;

		// compute distance between new and old front
		int nmax = 2;
		labelHashSet usedCells;

		for (labelHashSet::iterator i = newFront.begin(); i != newFront.end(); i++)
		{
			label cellI = i.key();

			// find nmax closest cells to cellI
			HashTable<scalar,label> cellDist(nmax);
			getNeighboursRecursive(mesh, cellI, nmax, 4, frozenCells, cellDist);

			if (cellDist.size() == 1)
			{
				// only one closest cell
				label cellI1 = cellDist.begin().key();
				d[cellI] = d[cellI1] + cellDist[cellI1];
				g[cellI] = (c[cellI] - c[cellI1])/cellDist[cellI1];
			}
			else if (cellDist.size() >= 2)
			{
				// two close cells
				// get the cell indices
				HashTable<scalar,label>::iterator k = cellDist.begin();
				label cellI1 = k.key(); k++;
				label cellI2 = k.key();

				// get center and gradients
				vector p = c[cellI];
				vector g1 = g[cellI1];
				vector g2 = g[cellI2];

				// get the points of equal mean distance
				scalar dmean = 0.5*(d[cellI1] + d[cellI2]);
				vector p1 = c[cellI1] + (dmean - d[cellI1])*g1;
				vector p2 = c[cellI2] + (dmean - d[cellI2])*g2;
				vector p12 = p1 - p2;
				scalar p12g1 = p12 & g1;
				scalar p12g2 = p12 & g2;
				scalar g12 = g1 & g2;

				// compute interpolation parameter t
				// i.e. min distance of p to
				// x = p1 + t*(p2-p1)
				// i.e. (x - p)*p12 = 0
				scalar t = ((p1 - p) & p12)/(p12 & p12);

				// discriminant g1 and g2 parallel
				scalar D = 1.0 - g12*g12;
				
				if (fabs(D) < SMALL)
				{
					// this case rarely happens
					Info << "parallel" << endl;

					if (g12 > 0) {
						// lines are parallel
						g[cellI] = g1*(1-t) + g2*t;
						g[cellI] /= mag(g[cellI]);
						vector q = p1*(1-t) + p2*t;
						d[cellI] = dmean + mag(p - q);
					}
					else {
						// lines are anti parallel
						scalar d12 = mag(p12);
						d[cellI] = (d[cellI1] + t*d12)*(1-t) + (d[cellI2] + (1-t)*d12)*t;
						g[cellI] = vector(0, 0, 0);
					}
				}
				else
				{
					// find the intersection/min. distance point of the lines
					// line1: x1 = p1 + t1*g1
					// line2: x2 = p2 + t2*g2
					// i.e. (x1-x2)*g1 = 0 and (x1-x2)*g2 = 0
					// i.e. p12g1 + t1 - t2*g12 = 0
					//      p12g2 + t1*g12 - t2 = 0

					// i.e. p12g1 + t1  -p12g2*g12 - t1*g12*g12 = 0
					// compute minimum distance points x1 and x2
					scalar t1 = (p12g2*g12 - p12g1)/D;
					scalar t2 = p12g2 + t1*g12;
					vector x1 = p1 + t1*g1;
					vector x2 = p2 + t2*g2;

					vector d1 = p - x1;
					vector d2 = p - x2;
					vector n1 = d1/mag(d1);
					vector n2 = d2/mag(d2);

					// NOTE: maybe better 0.5*(n1 + n2)?
					g[cellI] = n1*(1-t) + n2*t;
					g[cellI] /= mag(g[cellI]);

					if ((g[cellI] & (g1 + g2)) < 0) {
						g[cellI] *= -1.0;
					}
					
					vector q1 = x1 - n1*t1;
					vector q2 = x2 - n2*t2;

					// NOTE: maybe better 0.5*(q1 + q2)?
					vector q = q1*(1-t) + q2*t;

					d[cellI] = dmean + mag(p-q)*(((p - q) & g[cellI]) < 0 ? -1 : 1);
				}
			}

//			for (HashTable<scalar,label>::iterator k = cellDist.begin(); k != cellDist.end(); k++) {
//				usedCells.set(k.key());
//			}
		}

		// remove unused cells from frozen cell set
//		labelHashSet unUsedCells = frozenCells;
//		unUsedCells.erase(usedCells);
//		frozenCells.erase(unUsedCells);

		// add new front to frozen cell set
		for (labelHashSet::iterator i = newFront.begin(); i != newFront.end(); i++) {
			frozenCells.set(i.key());
		}

		Info << frozenCells.size() << endl;
	}

#else
	const vectorField& c =	mesh.cellCentres();

	labelHashSet frozenCells;
	labelHashSet newFront;

		frozenCells.insert(23415);

	// set inital distance value
	for (labelHashSet::iterator i = frozenCells.begin(); i != frozenCells.end(); i++) {
		d[i.key()] = 0.0;
	}

	for(;;)
	{
		// clear current front
		newFront.clear();

		for (labelHashSet::iterator i = frozenCells.begin(); i != frozenCells.end(); i++)
		{
			// get cell neighbours
			label cellI = i.key();
			const labelList& neighbours = mesh.cellCells()[cellI];

			// iterate over cell neighbours
			for (int j = 0; j < neighbours.size(); j++)
			{
				label neighbourCellI = neighbours[j];

				// skip already computed cells
				if (d[neighbourCellI] >= 0.0) {
					continue;
				}
				newFront.set(neighbourCellI);
			}
		}

		if (newFront.size() == 0) break;

		// compute distance between new and old front
		int nmax = 2;
		labelHashSet usedCells;

		for (labelHashSet::iterator i = newFront.begin(); i != newFront.end(); i++)
		{
			label cellI = i.key();

			// find nmax closest cells to cellI
			HashTable<scalar,label> cellDist(nmax);
			getNeighboursRecursive(mesh, cellI, nmax, 4, frozenCells, cellDist);

			if (cellDist.size() == 1)
			{
				// only one closest cell
				label cellI1 = cellDist.begin().key();
				d[cellI] = d[cellI1] + cellDist[cellI1];
			}
			else if (cellDist.size() >= 2)
			{
				// two close cells
				HashTable<scalar,label>::iterator k = cellDist.begin();
				label cellI1 = k.key(); k++;
				label cellI2 = k.key();
				vector d12 = c[cellI2] - c[cellI1];
				scalar alpha = ((c[cellI] - c[cellI1]) & d12)/(d12 & d12);
				vector dmin = (c[cellI] - c[cellI1]) - alpha*d12;
//				d[cellI] = (d[cellI1]*cellDist[cellI1] + d[cellI2]*cellDist[cellI2])/(cellDist[cellI1] + cellDist[cellI2]) + mag(dmin);
				d[cellI] = d[cellI1] + alpha*(d[cellI2] - d[cellI1]) + mag(dmin);
			}
			else if (cellDist.size() >= 3)
			{
/*
				HashTable<scalar,label>::iterator k = cellDist.begin();
				label cellI1 = k.key(); k++;
				label cellI2 = k.key(); k++;
				label cellI3 = k.key();
				vector d1x = c[cellI] - c[cellI1];
				vector d21 = c[cellI2] - c[cellI1];
				vector d31 = c[cellI1] - c[cellI1];
				vector n = (d21 ^ d31);
*/
			}

			for (HashTable<scalar,label>::iterator k = cellDist.begin(); k != cellDist.end(); k++) {
				usedCells.insert(k.key());
			}
		}

		// remove unused cells from frozen cell set
		labelHashSet unUsedCells = frozenCells;
		unUsedCells.erase(usedCells);
		frozenCells.erase(unUsedCells);

		// add new front to frozen cell set
		for (labelHashSet::iterator i = newFront.begin(); i != newFront.end(); i++) {
			frozenCells.set(i.key());
		}

		Info << frozenCells.size() << endl;
	}
#endif

	// write distance function
	d.write();

	// write gradient function
	g.write();

	Info<< "End\n" << endl;

	return 0;
}


// Get maximum nmax neighbour cells around cellI recursively (levels deep), which are also contained in the constraint set.
// Closer cells are preferred (i.e. nmax closest cells which are in the constraint set are returned).
// The result hash contains the cell ids and the center-center distance to cellI.
void getNeighboursRecursive(polyMesh& mesh, label cellI, int nmax, int levels, const labelHashSet& constraint, HashTable<scalar,label>& result, label currentCellI, labelHashSet& visitedCells)
{
	// stop recursion if necessary
	if (levels <= 0) return;

	// get list of neighbour cells for cellI
	const labelList& neighbours = mesh.cellCells()[currentCellI];

	// get vector field of cell centers
	const vectorField& center = mesh.cellCentres();

	// iterate over cell neighbours
	for (int j = 0; j < neighbours.size(); j++)
	{
		// get neighbour cell index
		label neighbourCellI = neighbours[j];

		// skip already visited cells
		if (visitedCells.found(neighbourCellI)) {
			continue;
		}

		// add cells only which are in the constraint set
		if (constraint.found(neighbourCellI))
		{
			// compute cell distance
			scalar dist = mag(center[cellI] - center[neighbourCellI]);

			if (result.size() < nmax) {
				// fill the result array up to nmax elements
				result.set(neighbourCellI, dist);
			}
			else
			{
				// find the element with maximum distance in results
				label kMax = result.begin().key();
				for (HashTable<scalar,label>::iterator k = result.begin(); k != result.end(); k++) {
					if (result[k.key()] > result[kMax]) {
						kMax = k.key();
					}
				}

				// replace the element with maximum distance by neighbourCellI, if closer
				if (dist < result[kMax]) {
					result.erase(kMax);
					result.set(neighbourCellI, dist);
				}
			}

		}

		// mark cell as visited
		visitedCells.set(neighbourCellI);

		// recursive call
		if (levels > 1) {
			getNeighboursRecursive(mesh, cellI, nmax, levels-1, constraint, result, neighbourCellI, visitedCells);
		}
	}
}

void getNeighboursRecursive(polyMesh& mesh, label cellI, int nmax, int levels, const labelHashSet& constraint, HashTable<scalar,label>& result)
{
	labelHashSet visitedCells(8*levels*levels);
	visitedCells.set(cellI);
	getNeighboursRecursive(mesh, cellI, nmax, levels, constraint, result, cellI, visitedCells);
}

// ************************************************************************* //


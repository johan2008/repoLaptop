#include "R3Shapes/R3Shapes.h"
#include "RNExt.h"
////////////////////////////////////////////////////////////////////////
// Structures
////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////
// Input functions
////////////////////////////////////////////////////////////////////////

// Mesh
R3Mesh *
ReadMesh(char *filename)
{
  	// Start statistics
  	RNTime start_time;
  	start_time.Read();

  	// Allocate mesh
  	R3Mesh *mesh = new R3Mesh();
  	assert(mesh);

  	// Read mesh from file
  	if (!mesh->ReadFile(filename)) {
    	fprintf(stderr, "Unable to read mesh from %s\n", filename);
    	return NULL;
  	}

  	// Print statistics
/*  	if (print_verbose) {*/
        //printf("Read mesh from %s ...\n", filename);
        //printf("  Time = %.2f seconds\n", start_time.Elapsed());
        //printf("  # Faces = %d\n", mesh->NFaces());
        //printf("  # Edges = %d\n", mesh->NEdges());
        //printf("  # Vertices = %d\n", mesh->NVertices());
        //fflush(stdout);
  	/*}*/

  	// Return success
  	return mesh;
}

// Properties
R3MeshPropertySet *
ReadProperties(R3Mesh *mesh, const char *filename)
{
  	// Start statistics
  	RNTime start_time;
  	start_time.Read();

  	// Allocate properties
  	R3MeshPropertySet *properties = new R3MeshPropertySet(mesh);
  	if (!properties) {
    	fprintf(stderr, "Unable to allocate properties for %s\n", filename);
    	return NULL;
  	}

  	// Read properties from file
  	if (!properties->Read(filename)) {
    	delete properties;
    	return NULL;
  	}

  	// Print statistics
/*  	if (print_verbose) {*/
        //printf("Read properties from %s ...\n", filename);
        //printf("  Time = %.2f seconds\n", start_time.Elapsed());
        //printf("  # Properties = %d\n", properties->NProperties());
        //printf("  # Vertices = %d\n", properties->Mesh()->NVertices());
        //fflush(stdout);
  	/*}*/

  	// Return property set
  	return properties;
}

// Property
R3MeshProperty *
ReadProperty(R3Mesh *mesh, const char *filename)
{
	  // Start statistics
  	RNTime start_time;
  	start_time.Read();

  	// Allocate properties
 	R3MeshProperty *property = new R3MeshProperty(mesh);
  	if (!property) {
    	fprintf(stderr, "Unable to allocate property for %s\n", filename);
    	return NULL;
  	}

  	// Read properties from file
  	if (!property->Read(filename)) {
    	delete property;
    	return NULL;
  	}

  	// Print statistics
/*  	if (print_verbose) {*/
        //printf("Read property from %s ...\n", filename);
        //printf("  Time = %.2f seconds\n", start_time.Elapsed());
        //printf("  # Vertices = %d\n", property->Mesh()->NVertices());
        //fflush(stdout);
  	/*}*/

  	// Return property set
  	return property;
}

// Points
RNArray<R3MeshVertex*>*
ReadPoints(R3Mesh* mesh, const char* filename)
{
	RNArray<R3MeshVertex*>* axis = new RNArray<R3MeshVertex*>;
	ifstream fin(filename);
	if (!fin)
	{
		cout<<"Cannot open point file "<<filename<<endl;
		exit(-1);
	}
	int temp;
	while (fin>>temp)
	{
		axis->Insert(mesh->Vertex(temp));
	}
	fin.close();
	return axis;
}

// Sparse Correspondence
SparseVertexCorrespondence *
ReadSparseVertexCorrespondence(R3Mesh *mesh0, R3Mesh *mesh1, char *filename)
{
  	// Start statistics
  	RNTime start_time;
  	start_time.Read();

  	// Parse filename extension
  	const char *extension;
  	if (!(extension = strrchr(filename, '.'))) {
    	printf("Filename %s has no extension (e.g., .cor)\n", filename);
    	return 0;
  	}

  	// Open file
  	FILE *fp = fopen(filename, "r");
  	if (!fp) {
    	fprintf(stderr, "Unable to open correspondence file %s\n", filename);
   	 	return NULL;
  	}

  	// Create correspondence
  	SparseVertexCorrespondence *correspondence = new SparseVertexCorrespondence(mesh0, mesh1);
  	if (!correspondence) {
    	fprintf(stderr, "Unable to allocate correspondence for %s\n", filename); 
    	return NULL; 
  	}

  	// Check filename extension
  	if (!strcmp(extension, ".cor")) {
    	// Read correspondences
    	int id0, id1;
    	while (fscanf(fp, "%d%d", &id0, &id1) == (unsigned int) 2) {
      		R3MeshVertex *vertex0 = mesh0->Vertex(id0);
      		R3MeshVertex *vertex1 = mesh1->Vertex(id1);
      		correspondence->vertices[0].Insert(vertex0);
      		correspondence->vertices[1].Insert(vertex1);
			correspondence->nvertices ++;
    	}
  	}
  	else if (!strcmp(extension, ".vk")) {
    	// Read header
    	double fdummy;
    	char format[256], dummy[256];
    	int nmeshes, ncorrespondences, idummy;
    	if (fscanf(fp, "%s%s%s%s%s%d%s%s%s%d%s%d%d%d%d%d", dummy, format,  dummy, dummy,  dummy, &nmeshes,
               dummy, dummy,  dummy, &idummy,   dummy, &idummy, &idummy, &idummy, &idummy, &ncorrespondences) != (unsigned int) 16) {
      		fprintf(stderr, "Unable to read %s\n", filename);
      		return NULL;
    	}

    	// Read correspondences
    	int id0, id1;
    	int count = 0;
    	while (fscanf(fp, "%d%lg%d%lg", &id0, &fdummy, &id1, &fdummy) == (unsigned int) 4) {
      		R3MeshVertex *vertex0 = mesh0->Vertex(id0);
      		R3MeshVertex *vertex1 = mesh1->Vertex(id1);
      		correspondence->vertices[0].Insert(vertex0);
      		correspondence->vertices[1].Insert(vertex1);
      		count++;
			correspondence->nvertices ++;
    	}

    	// Check number of correspondences
    	if (count != ncorrespondences) {
      	fprintf(stderr, "Mismatching number of correspondences in %s\n", filename);
      	return NULL;
    	}
  	}
  	else {
    	fprintf(stderr, "Unrecognized correspondence file extension: %s\n", extension);
    	return NULL;
  	}

  	// Close file
  	fclose(fp);

/*  	// Print statistics*/
      //if (print_verbose) {
        //printf("Read sparse correspondences from %s ...\n", filename);
        //printf("  Time = %.2f seconds\n", start_time.Elapsed());
        //printf("  # Correspondences = %d\n", correspondence->vertices[0].NEntries());
        //fflush(stdout);
  	/*}*/

  	// Return correspondence
  	return correspondence;
}

// Dense Correspondence
DenseVertexCorrespondence *
ReadDenseVertexCorrespondence(R3Mesh* mesh0, R3Mesh* mesh1, char* filename)
{
	// Create Dense Correspondence
	DenseVertexCorrespondence* map = new DenseVertexCorrespondence(mesh0, mesh1);

	// Open file
	ifstream fin(filename);
	if (!fin)
	{
		printf("Cannot open file for reading dense correspondence : %s\n", filename);
		delete map;
		return NULL;
	}
	int temp;
	int counter = 0;
	while (fin >> temp)
	{
		if (counter >= mesh0->NVertices())
		{
			printf("Map file does not match to Mesh : %s\n", filename);
			delete map;
			fin.close();
			return NULL;
		}
		map->vertices[counter] = mesh1->Vertex(temp);
		counter ++;
	}

	if (counter != mesh0->NVertices())
	{
		printf("Map file does not match Mesh : %s\n", filename);
		delete map;
		fin.close();
		return NULL;
	}
	
	fin.close();
	return map;
}

////////////////////////////////////////////////////////////////////////
// Output functions
////////////////////////////////////////////////////////////////////////

// Sparse Correspondence
int
WriteSparseVertexCorrespondence(SparseVertexCorrespondence *correspondence, char *filename)
{
	// Start statistics
  	RNTime start_time;
	start_time.Read();

	// Open file
	FILE *fp = fopen(filename, "w");
	if (!fp) {
		fprintf(stderr, "Unable to open sparse file %s\n", filename);
		return 0;
	}

	// Write correspondences
	for (int i=0; i < correspondence->nvertices; i++) {
		fprintf(fp, "%6d\t%6d\n", correspondence->mesh[0]->VertexID(correspondence->vertices[0].Kth(i)), 
				correspondence->mesh[1]->VertexID(correspondence->vertices[1].Kth(i)));
	}

	// Close file
	fclose(fp);

   /* // Print statistics*/
	//if (print_verbose) {
        //printf("Wrote sparse correspondences to %s ...\n", filename);
        //printf("  Time = %.2f seconds\n", start_time.Elapsed());
        //printf("  # Correspondences = %d\n", correspondence->nvertices);
        //fflush(stdout);
  	/*}*/

	// Return success
	return 1;	
}

// Dense Correspondence
int
WriteDenseVertexCorrespondence(DenseVertexCorrespondence *map, char *filename)
{
  	// Start statistics
  	RNTime start_time;
  	start_time.Read();

  	// Open file
  	FILE *fp = fopen(filename, "w");
  	if (!fp) {
    	fprintf(stderr, "Unable to open map file %s\n", filename);
    	return 0;
	}

	// Write correspondences
	for (int i = 0; i < map->mesh[0]->NVertices(); i++) {
		if (map->vertices[i] == NULL) fprintf(fp, "-1\n");
    	else fprintf(fp, "%d\n", map->mesh[1]->VertexID(map->vertices[i]));
  	}

  	// Close file
  	fclose(fp);

  	// Print statistics
/*  	if (print_verbose) {*/
        //printf("Wrote dense correspondences to %s ...\n", filename);
        //printf("  Time = %.2f seconds\n", start_time.Elapsed());
        //printf("  # Correspondences = %d\n", map->mesh[0]->NVertices());
        //fflush(stdout);
  	/*}*/

  	// Return success
  	return 1;
}

// Points
int
WritePoints(R3Mesh* mesh, RNArray<R3MeshVertex*> *points, const char* filename)
{
	ofstream fout(filename);
	if (!fout)
	{
		printf("Cannot open file for writing Points : %s\n", filename);
		return 0;
	}
	for (int i=0; i<points->NEntries(); i++)
		fout<<mesh->VertexID(points->Kth(i))<<endl;

	// Close file
	fout.close();
	return 1;
}


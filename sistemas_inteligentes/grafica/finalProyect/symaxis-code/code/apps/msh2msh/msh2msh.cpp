// Source file for the mesh converter program



// Include files 

#include "R3Shapes/R3Shapes.h"



// Program arguments

char *input_name = NULL;
char *output_name = NULL;
int flip_faces = 0;
int clean = 0;
RNScalar scale_factor = 0;
RNLength min_edge_length = 0;
RNLength max_edge_length = 0;
char *xform_name = NULL;
int scale_by_area = 0;
int scale_by_radius = 0;
int align_by_pca = 0;
int print_verbose = 0;



////////////////////////////////////////////////////////////////////////
// I/O STUFF
////////////////////////////////////////////////////////////////////////

static R3Mesh *
ReadMesh(char *mesh_name)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();

  // Allocate mesh
  R3Mesh *mesh = new R3Mesh();
  assert(mesh);

  // Read mesh from file
  if (!mesh->ReadFile(mesh_name)) {
    delete mesh;
    return NULL;
  }

  // Check if mesh is valid
  assert(mesh->IsValid());

  // Print statistics
  if (print_verbose) {
    printf("Read mesh ...\n");
    printf("  Time = %.2f seconds\n", start_time.Elapsed());
    printf("  # Faces = %d\n", mesh->NFaces());
    printf("  # Edges = %d\n", mesh->NEdges());
    printf("  # Vertices = %d\n", mesh->NVertices());
    fflush(stdout);
  }

  // Return success
  return mesh;
}



static int
WriteMesh(R3Mesh *mesh, char *filename)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();

  // Write mesh to file
  if (!mesh->WriteFile(filename)) return 0;

  // Print statistics
  if (print_verbose) {
    printf("Wrote mesh to %s ...\n", filename);
    printf("  Time = %.2f seconds\n", start_time.Elapsed());
    printf("  # Faces = %d\n", mesh->NFaces());
    printf("  # Edges = %d\n", mesh->NEdges());
    printf("  # Vertices = %d\n", mesh->NVertices());
    fflush(stdout);
  }

  // Return success
  return 1;
}



static int
ReadMatrix(R4Matrix& m, const char *filename)
{
  // Open file
  FILE *fp = fopen(filename, "r");
  if (!fp) {
    fprintf(stderr, "Unable to open matrix file: %s\n", filename);
    return 0;
  }

  // Read matrix from file
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      double value;
      fscanf(fp, "%lf", &value);
      m[i][j] = value;
    }
  }

  // Close file
  fclose(fp);

  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////
// PROGRAM ARGUMENT PARSING
////////////////////////////////////////////////////////////////////////

int ParseArgs(int argc, char **argv)
{
  // Check number of arguments
  if (argc == 1) {
    printf("Usage: mesh2mesh inputname outputname [options]\n");
    exit(0);
  }

  // Parse arguments
  argc--; argv++;
  while (argc > 0) {
    if ((*argv)[0] == '-') {
      if (!strcmp(*argv, "-v")) print_verbose = 1;
      else if (!strcmp(*argv, "-flip")) flip_faces = 1;
      else if (!strcmp(*argv, "-clean")) clean = 1;
      else if (!strcmp(*argv, "-align_by_pca")) align_by_pca = 1;
      else if (!strcmp(*argv, "-scale_by_area")) scale_by_area = 1;
	  else if (!strcmp(*argv, "-scale_by_radius")) scale_by_radius = 1;
      else if (!strcmp(*argv, "-scale")) { argv++; argc--; scale_factor = atof(*argv); }
      else if (!strcmp(*argv, "-xform")) { argv++; argc--; xform_name = *argv; }
      else if (!strcmp(*argv, "-min_edge_length")) { argv++; argc--; min_edge_length = atof(*argv); }
      else if (!strcmp(*argv, "-max_edge_length")) { argv++; argc--; max_edge_length = atof(*argv); }
      else { fprintf(stderr, "Invalid program argument: %s", *argv); exit(1); }
      argv++; argc--;
    }
    else {
      if (!input_name) input_name = *argv;
      else if (!output_name) output_name = *argv;
      else { fprintf(stderr, "Invalid program argument: %s", *argv); exit(1); }
      argv++; argc--;
    }
  }

  // Check input filename
  if (!input_name) {
    fprintf(stderr, "You did not specify an input file name.\n");
    return 0;
  }

  // Check output filename
  if (!output_name) {
    fprintf(stderr, "You did not specify an output file name.\n");
    return 0;
  }

  // Return OK status 
  return 1;
}



////////////////////////////////////////////////////////////////////////
// MAIN
////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
  // Check number of arguments
  if (!ParseArgs(argc, argv)) exit(1);

  // Read mesh
  R3Mesh *mesh = ReadMesh(input_name);
  if (!mesh) exit(-1);

  // Clean 
  if (clean) {
    mesh->DeleteUnusedEdges();
    mesh->DeleteUnusedVertices();
  }

  // Compute scale factor based on mesh area
  if ((scale_factor == 0) && (scale_by_area)) {
    RNArea area = mesh->Area();
    if (area > 0) scale_factor = 1 / sqrt(area);
  }

  //  Compute scale factor based on mesh radius
	if ((scale_factor == 0) && (scale_by_radius)) {
   		R3Point centroid = mesh->Centroid();
   		RNScalar r = 0;
   		for (int i=0; i<mesh->NVertices(); i++)
   		{
  			RNScalar d = R3Distance(mesh->VertexPosition(mesh->Vertex(i)), centroid);
 			if (d > r)
				r = d;
  		}
   		if (r > 0) scale_factor = 1 / r;
   	}
   
  // Normalize translation, rotation, and scale
  if (align_by_pca) {
    R3Affine xf = mesh->PCANormalizationTransformation();
    mesh->Transform(xf);
  }

  // Transform every vertex
  if (xform_name) {
    // Read xform
    R4Matrix m;
    if (ReadMatrix(m, xform_name)) {
      mesh->Transform(R3Affine(m));
    }
  }

  // Scale every vertex
  if ((scale_factor != 0) && (scale_factor != 1)) {
    for (int i = 0; i < mesh->NVertices(); i++) {
      R3MeshVertex *v = mesh->Vertex(i);
      R3Point p = mesh->VertexPosition(v);
      mesh->SetVertexPosition(v, p * scale_factor);
    }
  }

  // Flip every face
  if (flip_faces) {
    for (int i = 0; i < mesh->NFaces(); i++) {
      mesh->FlipFace(mesh->Face(i));
    }
  }

  // Subdivice edges that are too long
  if (max_edge_length > 0) {
    mesh->SubdivideLongEdges(max_edge_length);
  }

  // Split edges that are too long
  if (min_edge_length > 0) {
    mesh->CollapseShortEdges(min_edge_length);
  }

  // Write mesh
  if (!WriteMesh(mesh, output_name)) exit(-1);

  // Return success 
  return 0;
}


















// Source file 



// Include files 

#include "R3Shapes/R3Shapes.h"
#include <fstream>
#include <vector>
using namespace std;


////////////////////////////////////////////////////////////////////////
// Program arguments
////////////////////////////////////////////////////////////////////////

static char *input_mesh_name = NULL;
static char *input_points_name = NULL;
static char *output_properties_name = NULL;
static double start_radius[2] = {0.05, 0.05};
static double end_radius[2] = {0.25, 0.05};
static int nradii[2] = {2, 2};
static int norm_across = 0;
static int norm_along = 0;
static int blur_across = 0;
static int blur_along = 0;
static double sigma = 0.025;
static int print_verbose = 0;
static int print_debug = 0;

static double g_minimum = -FLT_MAX;
static double g_maximum = FLT_MAX;



////////////////////////////////////////////////////////////////////////
// Input/output
////////////////////////////////////////////////////////////////////////

static R3Mesh *
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
    fprintf(stderr, "Unable to read mesh %s\n", filename);
    return NULL;
  }

  // Normalize mesh scale by sqrt(area)
  RNScalar scale = sqrt(mesh->Area());
  for (int i = 0; i < mesh->NVertices(); i++) {
    R3MeshVertex *vertex = mesh->Vertex(i);
    R3Point position = mesh->VertexPosition(vertex);
    mesh->SetVertexPosition(vertex, position / scale);
  }

  // Print statistics
  if (print_verbose) {
    printf("Read mesh from %s ...\n", filename);
    printf("  Time = %.2f seconds\n", start_time.Elapsed());
    printf("  # Faces = %d\n", mesh->NFaces());
    printf("  # Edges = %d\n", mesh->NEdges());
    printf("  # Vertices = %d\n", mesh->NVertices());
    fflush(stdout);
  }

  // Return success
  return mesh;
}



static RNArray<R3MeshVertex *> *
ReadVertices(R3Mesh *mesh, const char *filename)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();

  // Allocate array of vertices
  RNArray<R3MeshVertex *> *vertices = new RNArray<R3MeshVertex *>();
  if (!vertices) {
    fprintf(stderr, "Unable to allocate array of vertices\n");
    return NULL;
  }

  // Open vertices file
  FILE *fp = fopen(filename, "r");
  if (!fp) {
    fprintf(stderr, "Unable to open vertices file: %s\n", filename);
    return NULL;
  }

  // Read vertices
  int vertex_index;
  while (fscanf(fp, "%d", &vertex_index) == (unsigned int) 1) {
    if ((vertex_index > 0) && (vertex_index < mesh->NVertices())) {
      R3MeshVertex *vertex = mesh->Vertex(vertex_index);
      vertices->Insert(vertex);
    }
  }

  // Close vertices file
  fclose(fp);

  // Print statistics
  if (print_verbose) {
    printf("Read point file: %s ...\n", filename);
    printf("  Time = %.2f seconds\n", start_time.Elapsed());
    printf("  # Vertices = %d\n", vertices->NEntries());
    fflush(stdout);
  }

  // Return vertices
  return vertices;
}



static int
WriteProperties(R3MeshPropertySet *properties, char *filename)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();

  // Write properties to file
  if (!properties->Write(filename)) {
    fprintf(stderr, "Unable to write properties to %s\n", filename);
    return 0;
  }

  // Print statistics
  if (print_verbose) {
    printf("Wrote properties to %s ...\n", filename);
    printf("  Time = %.2f seconds\n", start_time.Elapsed());
    printf("  # Properties = %d\n", properties->NProperties());
    printf("  # Vertices = %d\n", properties->Mesh()->NVertices());
    fflush(stdout);
  }

  // Return success
  return 1;
}



////////////////////////////////////////////////////////////////////////
// Axis creation
////////////////////////////////////////////////////////////////////////

static RNArray<R3MeshVertex *> *
CreateAxis(R3Mesh *mesh, const RNArray<R3MeshVertex *>& vertices)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();

  // Allocate array of vertices
  RNArray<R3MeshVertex *> *axis = new RNArray<R3MeshVertex *>();
  if (!axis) {
    fprintf(stderr, "Unable to allocate array of vertices\n");
    return NULL;
  }

  // Compute closed axis
  for (int i = 0; i < vertices.NEntries(); i++) {
    RNArray<R3MeshEdge *> shortest_axis_edges;
    R3MeshVertex *vertex0 = vertices.Kth(i);
    R3MeshVertex *vertex1 = vertices.Kth((i+1) % vertices.NEntries());
    mesh->DijkstraDistance(vertex1, vertex0, &shortest_axis_edges);
    R3MeshVertex *vertex = vertex0;
    for (int j = 0; j < shortest_axis_edges.NEntries(); j++) {
      R3MeshEdge *edge = shortest_axis_edges.Kth(j);
      assert(mesh->IsVertexOnEdge(vertex, edge));
      axis->Insert(vertex);
      vertex = mesh->VertexAcrossEdge(edge, vertex);
    }
  }

  // Print statistics
  if (print_verbose) {
    printf("Created axis ...\n");
    printf("  Time = %.2f seconds\n", start_time.Elapsed());
    printf("  # Vertices = %d\n", axis->NEntries());
    fflush(stdout);
  }

  // Return axis
  return axis;
}


////////////////////////////////////////////////////////////////////////
// Property calculation
////////////////////////////////////////////////////////////////////////

static R3MeshProperty *
CreateCurvaturesAlongAxis(R3Mesh *mesh, const RNArray<R3MeshVertex *>& axis, RNScalar radius)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();
  if (print_debug) {
    printf("  Creating curvatures along axis at radius %g ...\n", radius);
    fflush(stdout);
  }

  // Create property
  char buffer[1024];
  sprintf(buffer, "CurvatureAlongAxis_%g", radius);
  R3MeshProperty *property = new R3MeshProperty(mesh, buffer);
  if (!property) {
    fprintf(stderr, "Unable to allocate property\n");
    return NULL;
  }

  // Just checking
  if (axis.NEntries() < 3) return property;

  // Compute curvature at every vertex
  for (int k = 0; k < axis.NEntries(); k++) {
    R3MeshVertex *vertex = axis.Kth(k);
    R3Vector normal = mesh->VertexNormal(vertex);

    // Find first vertex
    R3MeshVertex *vertexA = axis[(k - 1 + axis.NEntries()) % axis.NEntries()];
    if (radius > 0) {
      R3MeshVertex *prev = vertexA;
      RNLength distance = R3Distance(mesh->VertexPosition(vertexA), mesh->VertexPosition(axis[k]));
      if (distance < radius) {
        for (int i = 2;  i < axis.NEntries()/3; i++) {
          vertexA = axis[(k - i + axis.NEntries()) % axis.NEntries()];
          distance += R3Distance(mesh->VertexPosition(vertexA), mesh->VertexPosition(prev));
          if (distance >= radius) break; 
          prev = vertexA;
        }
      }
    }

    // Find second vertex
    R3MeshVertex *vertexB = axis[(k + 1) % axis.NEntries()];
    if (radius > 0) {
      R3MeshVertex *prev = vertexB;
      RNLength distance = R3Distance(mesh->VertexPosition(vertexB), mesh->VertexPosition(axis[k]));
      if (distance < radius) {
        for (int i = 2;  i < axis.NEntries()/3; i++) {
          vertexB = axis[(k + i) % axis.NEntries()];
          distance += R3Distance(mesh->VertexPosition(vertexB), mesh->VertexPosition(prev));
          if (distance >= radius) break;
          prev = vertexB;
        }
      }
    }

    // Compute curvature
    R3Vector va = mesh->VertexPosition(vertexA) - mesh->VertexPosition(vertex);
    R3Vector vb = mesh->VertexPosition(vertexB) - mesh->VertexPosition(vertex);
    RNLength lena = va.Length();
    RNLength lenb = vb.Length();
    if (RNIsPositive(lena) && RNIsPositive(lenb)) {
      R3Vector vab = va + vb;
      RNScalar angle = R3InteriorAngle(va, vb);
      RNScalar sign = (vab.Dot(normal) > 0) ? -1 : 1;
      RNScalar curvature = sign * (RN_PI - angle) / (lena + lenb);
      property->SetVertexValue(vertex, curvature);
    }
  }

  // Print statistics
  if (print_debug) {
    printf("    Time = %.2f seconds\n", start_time.Elapsed());
    printf("    # Vertices = %d\n", axis.NEntries());
    fflush(stdout);
  }

  // Return property
  return property;
}



static R3MeshProperty *
CreateCurvaturesAcrossAxis(R3Mesh *mesh, const RNArray<R3MeshVertex *>& axis, RNScalar radius)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();
  if (print_debug) {
    printf("  Creating curvatures across axis at radius %g ...\n", radius);
    fflush(stdout);
  }

  // Create property
  char buffer[1024];
  sprintf(buffer, "CurvatureAcrossAxis_%g", radius);
  R3MeshProperty *property = new R3MeshProperty(mesh, buffer);
  if (!property) {
    fprintf(stderr, "Unable to allocate property\n");
    return NULL;
  }

  // Just checking
  if (axis.NEntries() < 3) return property;

  // Compute curvature at every vertex
  for (int k = 0; k < axis.NEntries(); k++) {
    R3MeshVertex *vertex = axis.Kth(k);
    R3Vector normal = mesh->VertexNormal(vertex);

    // Find first vertex
    R3MeshVertex *vertexA = axis[(k - 1 + axis.NEntries()) % axis.NEntries()];
    if (radius > 0) {
      R3MeshVertex *prev = vertexA;
      RNLength distance = R3Distance(mesh->VertexPosition(vertexA), mesh->VertexPosition(axis[k]));
      if (distance < radius) {
        for (int i = 2;  i < axis.NEntries()/3; i++) {
          vertexA = axis[(k - i + axis.NEntries()) % axis.NEntries()];
          distance += R3Distance(mesh->VertexPosition(vertexA), mesh->VertexPosition(prev));
          if (distance >= radius) break;
          prev = vertexA;
        }
      }
    }

    // Find second vertex
    R3MeshVertex *vertexB = axis[(k + 1) % axis.NEntries()];
    if (radius > 0) {
      R3MeshVertex *prev = vertexB;
      RNLength distance = R3Distance(mesh->VertexPosition(vertexB), mesh->VertexPosition(axis[k]));
      if (distance < radius) {
        for (int i = 2;  i < axis.NEntries()/3; i++) {
          vertexB = axis[(k + i) % axis.NEntries()];
          distance += R3Distance(mesh->VertexPosition(vertexB), mesh->VertexPosition(prev));
          if (distance >= radius) break;
          prev = vertexB;
        }
      }
    }

    // Find vertices by tracing path from vertex in directions perpendicular to axis
    R3Point positionA = mesh->VertexPosition(vertexA);
    R3Point positionB = mesh->VertexPosition(vertexB);
    R3Vector direction = (positionB - positionA) % normal;
    if (RNIsZero(direction.Length())) continue;
    direction.Normalize();
    mesh->TracePath(vertex, direction, radius, &positionA);
    mesh->TracePath(vertex, -direction, radius, &positionB);

    // Compute curvature
    R3Vector va = positionA - mesh->VertexPosition(vertex);
    R3Vector vb = positionB - mesh->VertexPosition(vertex);
    RNLength lena = va.Length();
    RNLength lenb = vb.Length();
    if (RNIsPositive(lena) && RNIsPositive(lenb)) {
      R3Vector vab = va + vb;
      RNScalar angle = R3InteriorAngle(va, vb);
      RNScalar sign = (vab.Dot(normal) > 0) ? -1 : 1;
      RNScalar curvature = sign * (RN_PI - angle) / (lena + lenb);
      property->SetVertexValue(vertex, curvature);
    }
  }

  // Print statistics
  if (print_debug) {
    printf("    Time = %.2f seconds\n", start_time.Elapsed());
    printf("    # Vertices = %d\n", axis.NEntries());
    fflush(stdout);
  }

  // Return property
  return property;
}



////////////////////////////////////////////////////////////////////////
// Property processing
////////////////////////////////////////////////////////////////////////

static int
BlurPropertyAlongAxis(R3MeshProperty *property, const RNArray<R3MeshVertex *>& axis, RNLength sigma)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();
  if (print_debug) {
    printf("  Blurring %s with sigma %g ...\n", property->Name(), sigma);
    fflush(stdout);
  }

  // Get convenient variables
  R3Mesh *mesh = property->Mesh();
  RNScalar denom = 2 * sigma * sigma;
    
  // Blur values
  R3MeshProperty original_property(*property);
  RNLength radius = 3 * sigma;
  for (int k = 0; k < axis.NEntries(); k++) {
    R3MeshVertex *vertex = axis.Kth(k);

    // Initalize
    RNScalar total_value = original_property.VertexValue(vertex);
    RNScalar total_weight = 1;

    // Blending in vertices above
    if (sigma > 0) {
      for (int dir = -1; dir <= 1; dir += 2) {
        RNLength distance = 0;
        R3MeshVertex *prev = vertex;
        for (int i = 1;  i < axis.NEntries()/3; i++) {
          R3MeshVertex *cur = axis[(k + dir*i + axis.NEntries()) % axis.NEntries()];
          distance += R3Distance(mesh->VertexPosition(cur), mesh->VertexPosition(prev));
          RNScalar value = original_property.VertexValue(cur);
          RNScalar weight = exp(-distance * distance / denom);
          total_value += weight * value;
          total_weight += weight;
          if (distance >= radius) break;
          prev = cur;
        }
      }
    }

    // Assign vertex value
    property->SetVertexValue(vertex, total_value / total_weight);
  }

  // Print statistics
  if (print_debug) {
    printf("    Time = %.2f seconds\n", start_time.Elapsed());
    printf("    # Vertices = %d\n", axis.NEntries());
    fflush(stdout);
  }

  // Return success
  return 1;
}



static int
NormalizePropertyAlongAxis(R3MeshProperty *property, const RNArray<R3MeshVertex *>& axis)
{
  // Just checking
  if (axis.NEntries() == 0) return 1;

  // Start statistics
  RNTime start_time;
  start_time.Read();
  if (print_debug) {
    printf("  Normalizing %s ...\n", property->Name());
    fflush(stdout);
  }

  // Compute sum
  RNScalar sum = 0;
  for (int k = 0; k < axis.NEntries(); k++) {
    R3MeshVertex *vertex = axis.Kth(k);
    sum += property->VertexValue(vertex);
  }

  // Compute mean
  RNScalar mean = sum / axis.NEntries();

  // Compute rss
  RNScalar rss = 0;
  for (int k = 0; k < axis.NEntries(); k++) {
    R3MeshVertex *vertex = axis.Kth(k);
    RNScalar delta = property->VertexValue(vertex) - mean;
    rss += delta * delta;
  }

  // Compute stddev
  RNScalar variance = rss / axis.NEntries();
  RNScalar stddev = sqrt(variance);

  // Normalize
  if (stddev > 0) {
    for (int k = 0; k < axis.NEntries(); k++) {
      R3MeshVertex *vertex = axis.Kth(k);
      RNScalar value = (property->VertexValue(vertex) - mean) / stddev;
      property->SetVertexValue(vertex, value);
    }
  }

  // Print statistics
  if (print_debug) {
    printf("    Time = %.2f seconds\n", start_time.Elapsed());
    printf("    # Vertices = %d\n", axis.NEntries());
    fflush(stdout);
  }

  // Return success
  return 1;
}



static R3MeshPropertySet *
CreateProperties(R3Mesh *mesh, const RNArray<R3MeshVertex *>& axis) 
{
  // Start statistics
  RNTime start_time;
  start_time.Read();
  if (print_verbose) {
    printf("Creating properties ...\n");
    fflush(stdout);
  }


  // Create property set
  R3MeshPropertySet *properties = new R3MeshPropertySet(mesh);
  if (!properties) {
    fprintf(stderr, "Unable to create properties\n");
    return NULL;
  }

  // Get convenient variables (along)
  RNLength step_radius = (nradii[0] > 1) ? (end_radius[0] - start_radius[0]) / (nradii[0] - 1) : 0;

  // Create curvatures along axis
  for (int i = 0; i < nradii[0]; i++) {
    RNLength radius = start_radius[0] + i * step_radius;
	if (print_verbose)
	{
	  printf("Radius: %0.3f ...\n", radius);
	}

    R3MeshProperty *property = CreateCurvaturesAlongAxis(mesh, axis, radius);
    properties->Insert(property);
  }

  step_radius = (nradii[1] > 1) ? (end_radius[1] - start_radius[1]) / (nradii[1] - 1) : 0;
  // Create curvatures across axis
  for (int i = 0; i < nradii[1]; i++) {
    RNLength radius = start_radius[1] + i * step_radius;
    R3MeshProperty *property = CreateCurvaturesAcrossAxis(mesh, axis, radius);
    properties->Insert(property);
  }

  // Blur properties
  for (int i = 0; i < properties->NProperties(); i++) {
    R3MeshProperty *property = properties->Property(i);
	// if it is along-axis
	if (i < nradii[0] && blur_along && sigma > 0)
	{
		int flag = BlurPropertyAlongAxis(property, axis, sigma);
		printf("Blur along-axis %d (sigma = %.3f) ...\n", i, sigma);
		if (!flag)
			return NULL;
	}
	else if (blur_across && sigma > 0)
	{
		printf("I guess it doesn't make sense to blur curvature across axis ...\n");
		exit(-1);
		int flag = BlurPropertyAlongAxis(property, axis, sigma);
		if (!flag)
			return NULL;
	}
  }

  // Normalize properties
  for (int i=0; i<properties->NProperties(); i++)
  {
	  R3MeshProperty * property = properties->Property(i);
	  // if it is along-axis
	  if (i < nradii[0] && norm_along)
	  {
		  printf("Normalize along-axis %d ...\n", i);
		  int flag = NormalizePropertyAlongAxis(property, axis);
		  if (!flag)
			  return NULL;
	  }
	  else if (norm_across)
	  {
		  printf("I guess it doesn't make sense to normalize curvature across axis ...\n");
//		  exit(-1);
		  int flag = NormalizePropertyAlongAxis(property, axis);
		  if (!flag)
			  return NULL;
	  }
  }
  // Print statistics
  if (print_verbose) {
    printf("  Time = %.2f seconds\n", start_time.Elapsed());
    printf("  # Properties = %d\n", properties->NProperties());
    fflush(stdout);
  }

  // Return properties
  return properties;
}



////////////////////////////////////////////////////////////////////////
// Parse program arguments
////////////////////////////////////////////////////////////////////////

static int 
ParseArgs(int argc, char **argv)
{
  // Parse arguments
  argc--; argv++;
  while (argc > 0) {
    if ((*argv)[0] == '-') {
      if (!strcmp(*argv, "-v")) print_verbose = 1;
      else if (!strcmp(*argv, "-debug")) print_debug = 1;
	  else if (!strcmp(*argv, "-norm_across"))
	  {
		  norm_across = 1;
	  }
	  else if (!strcmp(*argv, "-norm_along"))
	  {
		  norm_along = 1;
	  }
	  else if (!strcmp(*argv, "-blur_across"))
	  {
		  blur_across = 1;
	  }
	  else if (!strcmp(*argv, "-blur_along"))
	  {
		  blur_along = 1;
	  }
	  else if (!strcmp(*argv, "-clamp"))
	  {
		  argc--; argv++;
		  g_minimum = atof(*argv);
		  argc--; argv++;
		  g_maximum = atof(*argv);
	  }
//    else if (!strcmp(*argv, "-normalize")) normalize = 1;
	  else if (!strcmp(*argv, "-along"))
	  {
		  argc--; argv++;
		  nradii[0] = atoi(*argv);
		  argc--; argv++;
		  start_radius[0] = atof(*argv);
		  argc--; argv++;
		  end_radius[0] = atof(*argv);
	  }
	  else if (!strcmp(*argv, "-across"))
	  {
		  argc--; argv++;
		  nradii[1] = atoi(*argv);
		  argc--; argv++;
		  start_radius[1] = atof(*argv);
		  argc--; argv++;
		  end_radius[1] = atof(*argv);
	  }
   /*   else if (!strcmp(*argv, "-nradii"))*/
	  /*{*/
		  /*argc--; argv++;*/
		  /*nradii = atoi(*argv);*/
	  /*}*/
	  /*else if (!strcmp(*argv, "-start"))*/
	  /*{*/
		  /*argc--; argv++;*/
		  /*start_radius = atof(*argv);*/
	  /*}*/
	  /*else if (!strcmp(*argv, "-end"))*/
	  /*{*/
		  /*argc--; argv++;*/
		  /*end_radius = atof(*argv);*/
	  /*}*/
      else if (!strcmp(*argv, "-sigma")) { 
        argc--; argv++; sigma = atof(*argv); 
      }
/*      else if (!strcmp(*argv, "-radius")) { */
        /*argc--; argv++; start_radius = atof(*argv); */
        /*end_radius = start_radius;*/
        /*nradii = 1; */
      /*}*/
      else { 
        fprintf(stderr, "Invalid program argument: %s", *argv); 
        exit(1); 
      }
      argv++; argc--;
    }
    else {
      if (!input_mesh_name) input_mesh_name = *argv;
      else if (!input_points_name) input_points_name = *argv;
      else if (!output_properties_name) output_properties_name = *argv;
      else { fprintf(stderr, "Invalid program argument: %s", *argv); exit(1); }
      argv++; argc--;
    }
  }

  // Check input filename
  if (!input_mesh_name || !input_points_name || !output_properties_name) {
    fprintf(stderr, "Usage: axis2prp input_mesh input_points output_properties [-v, -debug,-sigma, -along n s e, -across n s e, -norm_across, -norm_along]\n");
    return 0;
  }

  // Return OK status 
  return 1;
}

static R3MeshVertex* 
NearestNeighbor(R3Mesh* mesh, R3Point pos, R3MeshSearchTree* search_tree)
{
	R3MeshIntersection closest;
	search_tree->FindClosest(pos, closest, 0, 1.0);
	if (closest.type == R3_MESH_NULL_TYPE)
	{
		return NULL;
	}

	RNLength closest_distance_to_surface = closest.t;
	R3Point closest_point_on_surface = closest.point;

	// Determine closest vertex
	R3MeshVertex *closest_vertex = NULL;
   	RNLength closest_distance_to_vertex = FLT_MAX;
   	if (closest.type == R3_MESH_VERTEX_TYPE) {
    	 // Closest point was on vertex
     	closest_vertex = closest.vertex;
     	closest_distance_to_vertex = closest.t;
   	}
   	else if (closest.type == R3_MESH_EDGE_TYPE) {
     	// Closest point was on edge
     	R3MeshVertex *vertex0 = mesh->VertexOnEdge(closest.edge, 0);
     	RNLength d0 = R3Distance(pos, mesh->VertexPosition(vertex0));
     	if (d0 < closest_distance_to_vertex) {
       		closest_vertex = vertex0;
       		closest_distance_to_vertex = d0;
     	}
     	R3MeshVertex *vertex1 = mesh->VertexOnEdge(closest.edge,1);
     	RNLength d1 = R3Distance(pos, mesh->VertexPosition(vertex1));
     	if (d1 < closest_distance_to_vertex) {
       		closest_vertex = vertex1;
       		closest_distance_to_vertex = d1;
     	}
   	}
   	else if (closest.type == R3_MESH_FACE_TYPE) {
     	// Closest point was in middle of face
     	R3MeshVertex *vertex0 = mesh->VertexOnFace(closest.face, 0);
     	RNLength d0 = R3Distance(pos, mesh->VertexPosition(vertex0));
     	if (d0 < closest_distance_to_vertex) {
       		closest_vertex = vertex0;
       		closest_distance_to_vertex = d0;
     	}
     	R3MeshVertex *vertex1 = mesh->VertexOnFace(closest.face,1);
     	RNLength d1 = R3Distance(pos, mesh->VertexPosition(vertex1));
     	if (d1 < closest_distance_to_vertex) {
       		closest_vertex = vertex1;
       		closest_distance_to_vertex = d1;
     	}
     	R3MeshVertex *vertex2 = mesh->VertexOnFace(closest.face,2);
     	RNLength d2 = R3Distance(pos, mesh->VertexPosition(vertex2));
     	if (d2 < closest_distance_to_vertex) {
       		closest_vertex = vertex2;
       		closest_distance_to_vertex = d2;
     	}
   	}

	return closest_vertex;
}


static vector<R3Point>*
ReadPositions(const char* filename)
{
	vector<R3Point>* points = new vector<R3Point>;
	ifstream fin(filename);
	RNScalar temp[3];
	while (fin>>temp[0]>>temp[1]>>temp[2])
	{
		R3Point p(temp[0], temp[1],temp[2]);
		points->push_back(p);
	}
	fin.close();
	return points;
}
static int
WriteAxisProperties(RNArray<R3MeshVertex*>* vertices, R3MeshPropertySet* properties, 
	 const char* filename)
{
//	R3Mesh* mesh = properties->Mesh();
	ofstream fout(filename);
	fout<<properties->NProperties()<<" "<<vertices->NEntries()<<endl;
	for (int k=0; k<properties->NProperties(); k++)
	{
		R3MeshProperty* prp = properties->Property(k);
		for (int i=0; i<vertices->NEntries(); i++)
		{
			R3MeshVertex* v = vertices->Kth(i);
			RNScalar val = prp->VertexValue(v);
			if (val < g_minimum)
			{
				val = g_minimum;
			}
			else if (val > g_maximum)
			{
				val = g_maximum;
			}
			fout<<val<<" ";
		}
		fout<<endl;
	}
	fout.close();
	return 1;
}

////////////////////////////////////////////////////////////////////////
// Main
////////////////////////////////////////////////////////////////////////

int 
main(int argc, char **argv)
{
  // Check number of arguments
  if (!ParseArgs(argc, argv)) exit(1);

  // Read mesh 
  R3Mesh *mesh = ReadMesh(input_mesh_name);
  if (!mesh) exit(-1);
  R3MeshSearchTree* search_tree = new R3MeshSearchTree(mesh);

  RNArray<R3MeshVertex *> *vertices = NULL;
  // Read vertices
  const char * extension = strrchr(input_points_name, '.');
  if (!strcmp(extension, ".pid"))
  {
	  printf("Reading *.pid points ...");
	vertices = ReadVertices(mesh, input_points_name);
  }
  else if (!strcmp(extension, ".xyz"))
  {
	  printf("Reading *.xyz points ...");
	 vector<R3Point>* points = ReadPositions(input_points_name);
	 vertices = new RNArray<R3MeshVertex*>;
	 for (int i=0; i<points->size(); i++)
	 {
		 R3MeshVertex* v = NearestNeighbor(mesh, (*points)[i], search_tree);
		 vertices->Insert(v);
	 }
  }
  if (!vertices) exit(-1);

  // Create path connecting vertices
  RNArray<R3MeshVertex *> *axis = CreateAxis(mesh, *vertices);
  if (!axis) exit(-1);

  // Create properties
  R3MeshPropertySet *properties = CreateProperties(mesh, *axis);
  if (!properties) exit(-1);

  // Write properties
  const char* _extension;
  if (!(_extension = strrchr(output_properties_name, '.'))) 
  {
	  printf("%s has no extension (e.g. parff or arff)", output_properties_name);
	  return 0;
  }
  if (!strcmp(_extension, ".arff"))
  {
	WriteProperties(properties, output_properties_name);

  }
  else if (!strcmp(_extension, ".parff"))
  {
	WriteAxisProperties(vertices, properties, output_properties_name);
  }

  // Return success 
  return 0;
}





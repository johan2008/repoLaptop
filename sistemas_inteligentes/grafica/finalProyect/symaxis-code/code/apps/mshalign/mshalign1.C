// Source file for the mesh alignment program



// Include files 

#include "R3Shapes/R3Shapes.h"



// Program arguments

static char *input1_name = NULL;
static char *input2_name = NULL;
static char *output1_name = NULL;
static int pca_translation = 1;
static int pca_scale = 1;
static int pca_rotation = 1; // 2 = only 180s
static int icp_translation = 1;
static int icp_scale = 1;
static int icp_rotation = 1;
static int max_points = 256;
static int print_verbose = 0;
static int print_debug = 0;



static int
CreatePoints(R3Mesh *mesh, R3Point *points, int max_points)
{
  // Check maximum number of points
  if (max_points <= 0) return 0;

  // Count total area of faces
  RNArea total_area = 0.0;
  for (int i = 0; i < mesh->NFaces(); i++) {
    R3MeshFace *face = mesh->Face(i);
    RNScalar face_area = mesh->FaceArea(face);
    mesh->SetFaceValue(face, face_area);
    total_area += face_area;
  }
    
  // Generate points with a uniform distribution over surface area
  int npoints = 0;
  RNSeedRandomScalar();
  for (int i = 0; i < mesh->NFaces(); i++) {
    R3MeshFace *face = mesh->Face(i);

    // Get vertex positions
    R3MeshVertex *v0 = mesh->VertexOnFace(face, 0);
    R3MeshVertex *v1 = mesh->VertexOnFace(face, 1);
    R3MeshVertex *v2 = mesh->VertexOnFace(face, 2);
    const R3Point& p0 = mesh->VertexPosition(v0);
    const R3Point& p1 = mesh->VertexPosition(v1);
    const R3Point& p2 = mesh->VertexPosition(v2);

    // Determine number of points for face 
    RNScalar ideal_face_npoints = max_points * mesh->FaceValue(face) / total_area;
    int face_npoints = (int) ideal_face_npoints;
    RNScalar remainder = ideal_face_npoints - face_npoints;
    if (remainder > RNRandomScalar()) face_npoints++;

    // Generate random points in face
    for (int j = 0; j < face_npoints; j++) {
      RNScalar r1 = sqrt(RNRandomScalar());
      RNScalar r2 = RNRandomScalar();
      RNScalar t0 = (1.0 - r1);
      RNScalar t1 = r1 * (1.0 - r2);
      RNScalar t2 = r1 * r2;
      points[npoints++] = t0*p0 + t1*p1 + t2*p2;
      if (npoints >= max_points) break;
    }

    // Check number of points created
    if (npoints >= max_points) break;
  }

#if 0
  if (0) {
    static int file_count = 1;
    char filename[256];
    sprintf(filename, "points%d.pts", file_count++);
    FILE *fp = fopen(filename, "wb");
    for (int i = 0; i < npoints; i++) {
      R3Point position = points[i];
      static float coordinates[6] = { 0 };
      coordinates[0] = points[i].X();
      coordinates[1] = points[i].Y();
      coordinates[2] = points[i].Z();
      fwrite(coordinates, sizeof(float), 6, fp);
    }
    fclose(fp);
  }
#endif

  // Return number of points actually created
  return npoints;
}



static int
CreatePointCorrespondences(
  R3Mesh *mesh1, const R3Point *points1, int npoints1,
  R3Mesh *mesh2, const R3Point *points2, int npoints2, 
  const R3Affine& affine12, const R3Affine& affine21,
  R3Point *correspondences1, R3Point *correspondences2, 
  int max_correspondences)
{
  // Initialize number of correspondences
  assert(max_correspondences == npoints1 + npoints2);
  int ncorrespondences = 0;

  // Compute correspondences for points1 -> mesh2
  static R3MeshSearchTree *tree2 = NULL;
  if (!tree2) tree2 = new R3MeshSearchTree(mesh2);
  for (int i = 0; i < npoints1; i++) {
    R3Point position1 = points1[i];
    position1.Transform(affine12);
    R3MeshIntersection closest;
    tree2->FindClosest(position1, closest);
    correspondences1[ncorrespondences] = points1[i];
    correspondences2[ncorrespondences] = closest.point;
    ncorrespondences++;
  }

  // Compute correspondences for points2 -> mesh1
  static R3MeshSearchTree *tree1 = NULL;
  if (!tree1) tree1 = new R3MeshSearchTree(mesh1);
  for (int i = 0; i < npoints2; i++) {
    R3Point position2 = points2[i];
    position2.Transform(affine21);
    R3MeshIntersection closest;
    tree1->FindClosest(position2, closest);
    correspondences1[ncorrespondences] = closest.point;
    correspondences2[ncorrespondences] = points2[i];
    ncorrespondences++;
  }

  // Return number of correspondences
  assert(ncorrespondences == npoints1 + npoints2);
  return ncorrespondences;
}



static RNScalar 
RMSD(R3Mesh *mesh1, const R3Point *points1, int npoints1,
     R3Mesh *mesh2, const R3Point *points2, int npoints2, 
     const R3Affine& affine12, const R3Affine& affine21)
{
  // Check number of points
  if (npoints1 == 0) return RN_INFINITY;
  if (npoints2 == 0) return RN_INFINITY;

  // Create array of correspondences
  int max_correspondences = npoints1 + npoints2;
  R3Point *correspondences1 = new R3Point [ max_correspondences ];
  R3Point *correspondences2 = new R3Point [ max_correspondences ];
  int ncorrespondences = CreatePointCorrespondences(
    mesh1, points1, npoints1, 
    mesh2, points2, npoints2, 
    affine12, affine21,
    correspondences1, correspondences2, 
    max_correspondences);

  // Add SSD of correspondences
  RNScalar ssd = 0;
  for (int i = 0; i < ncorrespondences; i++) {
    R3Point position1 = correspondences1[i];
    position1.Transform(affine12);
    R3Point& position2 = correspondences2[i];
    RNScalar d = R3Distance(position1, position2);
    ssd += d * d;
  }

  // Compute the RMSD
  RNScalar rmsd = sqrt(ssd / ncorrespondences);

  // Delete arrays of correspondences
  delete [] correspondences1;
  delete [] correspondences2;

  // Return RMSD
  return rmsd;
}



static int
ICPAlignmentTransformation(
  R3Mesh *mesh1, const R3Point *points1, int npoints1,
  R3Mesh *mesh2, const R3Point *points2, int npoints2, 
  R3Affine& affine12, R3Affine& affine21,
  int translation, int rotation, int scale)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();

  // Allocate arrays of correspondences
  int ncorrespondences = npoints1 + npoints2;
  R3Point *correspondences_buffer1 = new R3Point [2 * ncorrespondences];
  R3Point *correspondences_buffer2 = new R3Point [2 * ncorrespondences];
  assert(correspondences_buffer1 && correspondences_buffer2);

  // Iterate until max_iterations or converged
  RNBoolean converged = FALSE;
  const int max_iterations = 128;
  for (int iteration = 0; iteration < max_iterations; iteration++) {
    // Get arrays of correspondences
    R3Point *correspondences1 = &correspondences_buffer1[iteration%2 * ncorrespondences];
    R3Point *correspondences2 = &correspondences_buffer2[iteration%2 * ncorrespondences];

    // Update correspondences for aligning transformation
    CreatePointCorrespondences(
      mesh1, points1, npoints1, 
      mesh2, points2, npoints2, 
      affine12, affine21, 
      correspondences1, correspondences2, 
      ncorrespondences);

    // Update aligning transformation for correspondences
    R4Matrix matrix = R3AlignPoints(ncorrespondences, correspondences2, correspondences1, 
      NULL, translation, rotation, scale);
    affine12.Reset(matrix);
    affine21 = affine12.Inverse();

    // Check for convergence
    converged = FALSE;
    if (iteration > 0) {
      converged = TRUE;
      R3Point *prev_correspondences1 = &correspondences_buffer1[(1-(iteration%2)) * ncorrespondences];
      R3Point *prev_correspondences2 = &correspondences_buffer2[(1-(iteration%2)) * ncorrespondences];
      for (int i = 0; i < ncorrespondences; i++) {
        if (!R3Contains(correspondences1[i], prev_correspondences1[i])) { converged = FALSE; break; }
        if (!R3Contains(correspondences2[i], prev_correspondences2[i])) { converged = FALSE; break; }
      }
    }
    if (converged) break;
  }

  // Return whether converged
  return converged;
}



static int
Align(R3Mesh *mesh1, R3Mesh *mesh2, 
  int pca_translation, int pca_rotation, int pca_scale, 
  int icp_translation, int icp_rotation, int icp_scale)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();

  // Initialize transformation
  R3Affine affine12(R3identity_affine);
  R3Affine affine21(R3identity_affine);
  RNBoolean converged = TRUE;
  RNScalar rmsd = FLT_MAX;

  // Compute info for mesh1
  R3Point centroid1 = mesh1->Centroid();
  RNScalar radius1 = (pca_scale) ? mesh1->AverageRadius(&centroid1) : 1;
  RNScalar scale1 = (pca_scale && (radius1 > 0)) ? 1.0 / radius1 : 1;
  R3Triad axes1 = (pca_rotation==1) ? mesh1->PrincipleAxes(&centroid1) : R3xyz_triad;
  R3Point *points1 = new R3Point [ max_points ];
  int npoints1 = CreatePoints(mesh1, points1, max_points);
  if ((scale1 == 0) || (npoints1 == 0)) {
    fprintf(stderr, "Unable to process first mesh\n");
    return 0;
  }

  // Compute info for mesh2
  R3Point centroid2 = mesh2->Centroid();
  RNScalar scale2 = (pca_scale) ? mesh2->AverageRadius(&centroid2) : 1;
  R3Triad axes2 = (pca_rotation == 1) ? mesh2->PrincipleAxes(&centroid2) : R3xyz_triad;
  R3Affine affine02 = R3identity_affine;
  if (pca_translation) affine02.Translate(centroid2.Vector());
  if (pca_rotation==1) affine02.Transform(axes2.Matrix());
  if (pca_scale) affine02.Scale(scale2);
  R3Point *points2 = new R3Point [ max_points ];
  int npoints2 = CreatePoints(mesh2, points2, max_points);
  if ((scale2 == 0) || (npoints2 == 0)) {
    fprintf(stderr, "Unable to process first mesh\n");
    return 0;
  }

  // Compute RMSD for alignment with all flips of principle axes
  if (pca_rotation) {
    for (int dim1 = 0; dim1 < 3; dim1++) {
      for (int dim2 = 0; dim2 < 3; dim2++) {
        if (dim1 == dim2) continue;
        for (int dir1 = 0; dir1 < 2; dir1++) {
          for (int dir2 = 0; dir2 < 2; dir2++) {
            // Create triad of axes for flip
            R3Vector axis1a = (dir1 == 1) ? axes1[dim1] : -axes1[dim1];
            R3Vector axis1b = (dir2 == 1) ? axes1[dim2] : -axes1[dim2];
            R3Vector axis1c = axis1a % axis1b;
            R3Triad triad1(axis1a, axis1b, axis1c);

            // Compute transformation for mesh1 (with flip)
            R3Affine affine10 = R3identity_affine;
            if (pca_scale) affine10.Scale(scale1);
            affine10.Transform(triad1.InverseMatrix());
            if (pca_translation) affine10.Translate(-centroid1.Vector());

            // Compute composite transformation and its inverse
            R3Affine flipped_affine12 = affine02;
            flipped_affine12.Transform(affine10);
            R3Affine flipped_affine21 = flipped_affine12.Inverse();

            // Refine alignment with ICP
            RNBoolean flipped_converged = FALSE;
            if (icp_translation || icp_rotation || icp_scale) {
              flipped_converged = ICPAlignmentTransformation(
                mesh1, points1, npoints1, 
                mesh2, points2, npoints2, 
                flipped_affine12, flipped_affine21, 
                icp_translation, icp_rotation, icp_scale);
            }

            // Compute RMSD for transformation
            RNScalar flipped_rmsd = RMSD(
              mesh1, points1, npoints1, 
              mesh2, points2, npoints2, 
              flipped_affine12, flipped_affine21);

            // Check if best so far -- if so, save
            if (flipped_rmsd < rmsd) {
              affine12 = flipped_affine12;
              affine21 = flipped_affine21;
              converged = flipped_converged;
              rmsd = flipped_rmsd;
            }

            // Print statistics
            if (print_debug) {
              const R4Matrix& m = flipped_affine12.Matrix();
              printf("Computed alignment transformation for flip %d %d %d %d ...\n", dim1, dim2, dir1, dir2);
              printf("  Time = %.2f seconds\n", start_time.Elapsed());
              printf("  Max Points = %d\n", max_points);
              printf("  Matrix[0][0-3] = %g %g %g %g\n", m[0][0], m[0][1], m[0][2], m[0][3]);
              printf("  Matrix[1][0-3] = %g %g %g %g\n", m[1][0], m[1][1], m[1][2], m[1][3]);
              printf("  Matrix[2][0-3] = %g %g %g %g\n", m[2][0], m[2][1], m[2][2], m[2][3]);
              printf("  Matrix[3][0-3] = %g %g %g %g\n", m[3][0], m[3][1], m[3][2], m[3][3]);
              printf("  Scale = %g\n", flipped_affine12.ScaleFactor());
              printf("  Converged = %d\n", flipped_converged);
              printf("  RMSD = %g\n", flipped_rmsd);
              fflush(stdout);
            }
          }
        }
      }
    }
  }
  else {
    // Compute transformation for mesh1
    R3Affine affine10 = R3identity_affine;
    if (pca_scale) affine10.Scale(scale1);
    if (pca_translation) affine10.Translate(-centroid1.Vector());

    // Compute composite transformation
    affine12 = affine02;
    affine12.Transform(affine10);
    affine21 = affine12.Inverse();

    // Refine alignment with ICP
    if (icp_translation || icp_rotation || icp_scale) {
      converged = ICPAlignmentTransformation(
        mesh1, points1, npoints1, 
        mesh2, points2, npoints2, 
        affine12, affine21,
        icp_translation, icp_rotation, icp_scale);
    }

    // Compute RMSD
    rmsd = RMSD(
      mesh1, points1, npoints1, 
      mesh2, points2, npoints2, 
      affine12, affine21);
  }

  // Apply transformation
  mesh1->Transform(affine12);

  // Print statistics
  if (print_verbose) {
    const R4Matrix& m = affine12.Matrix();
    printf("Computed alignment transformation ...\n");
    printf("  Time = %.2f seconds\n", start_time.Elapsed());
    printf("  Matrix[0][0-3] = %g %g %g %g\n", m[0][0], m[0][1], m[0][2], m[0][3]);
    printf("  Matrix[1][0-3] = %g %g %g %g\n", m[1][0], m[1][1], m[1][2], m[1][3]);
    printf("  Matrix[2][0-3] = %g %g %g %g\n", m[2][0], m[2][1], m[2][2], m[2][3]);
    printf("  Matrix[3][0-3] = %g %g %g %g\n", m[3][0], m[3][1], m[3][2], m[3][3]);
    printf("  Scale = %g\n", affine12.ScaleFactor());
    printf("  Converged = %d\n", converged);
    printf("  RMSD = %g\n", rmsd);
    fflush(stdout);
  }

  // Return success
  return 1;
}

static int
Align(R3Mesh *mesh1, int translation, int rotation, int scale)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();

  // Compute transform to align mesh1 to canonical coordinate system
  R3Affine affine = mesh1->PCANormalizationTransformation(translation, (rotation==1), scale);

  // Apply transformation
  mesh1->Transform(affine);

  // Print statistics
  if (print_verbose) {
    const R4Matrix& m = affine.Matrix();
    printf("Computed alignment transformation with PCA ...\n");
    printf("  Time = %.2f seconds\n", start_time.Elapsed());
    printf("  Matrix[0][0-3] = %g %g %g %g\n", m[0][0], m[0][1], m[0][2], m[0][3]);
    printf("  Matrix[1][0-3] = %g %g %g %g\n", m[1][0], m[1][1], m[1][2], m[1][3]);
    printf("  Matrix[2][0-3] = %g %g %g %g\n", m[2][0], m[2][1], m[2][2], m[2][3]);
    printf("  Matrix[3][0-3] = %g %g %g %g\n", m[3][0], m[3][1], m[3][2], m[3][3]);
    printf("  Scale = %g\n", affine.ScaleFactor());
    fflush(stdout);
  }

  // Return success
  return 1;
}



static R3Mesh *
ReadMesh(const char *filename)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();

  // Allocate mesh
  R3Mesh *mesh = new R3Mesh();
  assert(mesh);

  // Read mesh from file
  if (!mesh->ReadFile(filename)) {
    delete mesh;
    return NULL;
  }

  // Check if mesh is valid
  assert(mesh->IsValid());

  // Print statistics
  if (print_verbose) {
    printf("Read mesh from %s ...\n", filename);
    printf("  Time = %.2f seconds\n", start_time.Elapsed());
    printf("  # Faces = %d\n", mesh->NFaces());
    printf("  # Edges = %d\n", mesh->NEdges());
    printf("  # Vertices = %d\n", mesh->NVertices());
    fflush(stdout);
  }

  // Return mesh
  return mesh;
}


static int
WriteMesh(R3Mesh *mesh, const char *filename)
{
  // Start statistics
  RNTime start_time;
  start_time.Read();

  // Write mesh to file
  if (!mesh->WriteFile(filename)) {
    return 0;
  }

  // Print statistics
  if (print_verbose) {
    printf("Wrote mesh to %s ...\n", filename);
    printf("  Time = %.2f seconds\n", start_time.Elapsed());
    printf("  # Faces = %d\n", mesh->NFaces());
    printf("  # Edges = %d\n", mesh->NEdges());
    printf("  # Vertices = %d\n", mesh->NVertices());
    fflush(stdout);
  }

  // Return succes
  return 1;
}


static int 
ParseArgs(int argc, char **argv)
{
  // Check number of arguments
  if (argc == 1) {
    printf("Usage: mshalign input1 [input2] output1 [options]\n");
    exit(0);
  }

  // Parse arguments
  argc--; argv++;
  while (argc > 0) {
    if ((*argv)[0] == '-') {
      if (!strcmp(*argv, "-v")) print_verbose = 1;
      else if (!strcmp(*argv, "-debug")) print_debug = 1;
      else if (!strcmp(*argv, "-no_pca")) { pca_translation = pca_rotation = pca_scale = 0; }
      else if (!strcmp(*argv, "-no_icp")) { icp_translation = icp_rotation = icp_scale = 0; }
      else if (!strcmp(*argv, "-no_translation")) { pca_translation = 0; icp_translation = 0; }
      else if (!strcmp(*argv, "-no_rotation")) { pca_rotation = 0; icp_rotation = 0; }
      else if (!strcmp(*argv, "-no_scale")) { pca_scale = 0; icp_scale = 0; }
      else if (!strcmp(*argv, "-axial_rotation")) { pca_rotation = 2; icp_rotation = 0; }
      else if (!strcmp(*argv, "-max_points")) { argc--; argv++; max_points = atoi(*argv); }
      else { fprintf(stderr, "Invalid program argument: %s", *argv); exit(1); }
      argv++; argc--;
    }
    else {
      if (!input1_name) input1_name = *argv;
      else if (!input2_name) input2_name = *argv;
      else if (!output1_name) output1_name = *argv;
      else { fprintf(stderr, "Invalid program argument: %s", *argv); exit(1); }
      argv++; argc--;
    }
  }

  // Fix output filename
  if (input2_name && !output1_name) {
    output1_name = input2_name;
    input2_name = NULL;
  }

  // Check input filename
  if (!input1_name || !output1_name) {
    printf("Usage: mshalign input1 [input2] output1 [options]\n");
    return 0;
  }

  // Return OK status 
  return 1;
}



int main(int argc, char **argv)
{
  // Check number of arguments
  if (!ParseArgs(argc, argv)) exit(1);

  // Read mesh1
  R3Mesh *mesh1 = ReadMesh(input1_name);
  if (!mesh1) exit(-1);

  // Check if aligning to second model
  if (input2_name) {
    // Read mesh2
    R3Mesh *mesh2 = ReadMesh(input2_name);
    if (!mesh2) exit(-1);

    // Align mesh1 to mesh2
    if (!Align(mesh1, mesh2, pca_translation, pca_rotation, pca_scale, icp_translation, icp_rotation, icp_scale)) exit(-1);
  }
  else {
    // Align mesh1 to canonical coordinate system
    if (!Align(mesh1, pca_translation, pca_rotation, pca_scale)) exit(-1);
  }

  // Output mesh1
  if (!WriteMesh(mesh1, output1_name)) exit(-1);

  // Return success 
  return 0;
}


















// Source file for the mesh property class



////////////////////////////////////////////////////////////////////////
// Include files 
////////////////////////////////////////////////////////////////////////

#include "R3Shapes/R3Shapes.h"



// Private variables

static int next_property_name = 1;



////////////////////////////////////////////////////////////////////////
// Constructors/destructors
////////////////////////////////////////////////////////////////////////

R3MeshProperty::
R3MeshProperty(R3Mesh *mesh, const char *name, const RNScalar *vertex_values)
  : mesh(mesh),
    nvalues(0),
    values(NULL),
    mean(RN_UNKNOWN),
    stddev(RN_UNKNOWN),
    minimum(RN_UNKNOWN),
    maximum(RN_UNKNOWN),
    median(RN_UNKNOWN),
    l2norm(RN_UNKNOWN)
{
  // Assign property name
  if (name) strncpy(this->name, name, 1024); 
  else sprintf(this->name, "Property%d", next_property_name++);
  this->name[1023] = '\0';
  
  // Assign vertex values
  this->nvalues = mesh->NVertices();
  this->values = new RNScalar [ this->nvalues ];
  if (vertex_values) {
    for (int i = 0; i < mesh->NVertices(); i++) {
      this->values[i] = vertex_values[i];
    }
  }
  else {
    for (int i = 0; i < mesh->NVertices(); i++) {
      this->values[i] = 0;
    }
  }
}



R3MeshProperty::
R3MeshProperty(const R3MeshProperty& property)
  : mesh(property.mesh),
    nvalues(0),
    values(NULL),
    mean(property.mean),
    stddev(property.stddev),
    minimum(property.minimum),
    maximum(property.maximum),
    median(property.median),
    l2norm(property.l2norm)
{
  // Copy property name
  strcpy(this->name, property.name);
  
  // Copy vertex values
  this->nvalues = mesh->NVertices();
  this->values = new RNScalar [ this->nvalues ];
  for (int i = 0; i < mesh->NVertices(); i++) {
    this->values[i] = property.values[i];
  }
}



R3MeshProperty::
~R3MeshProperty(void)
{
  // Delete vertex values
  if (values) delete [] values;
}



////////////////////////////////////////////////////////////////////////
// Statistics functions
////////////////////////////////////////////////////////////////////////

RNScalar R3MeshProperty::
Mean(void) const
{
  // Compute mean of property values
  if (mean == RN_UNKNOWN) {
    R3MeshProperty *property = (R3MeshProperty *) this;
    property->mean = 0;
    if (nvalues > 0) {
      // Compute sum 
      int count = 0;
      RNScalar sum = 0;
      for (int i = 0; i < nvalues; i++) {
        if (values[i] == RN_UNKNOWN) continue;
        sum += values[i];
        count++;
      }

      // Compute mean
      if (count > 0) property->mean = sum / count;
    }
  }

  // Return mean
  return mean;
}



RNScalar R3MeshProperty::
Minimum(void) const
{
  // Compute minimum of property values
  if (minimum == RN_UNKNOWN) {
    R3MeshProperty *property = (R3MeshProperty *) this;
    property->minimum = 0;
    if (nvalues > 0) {
      property->minimum = FLT_MAX;
      for (int i = 0; i < nvalues; i++) {
        if (values[i] == RN_UNKNOWN) continue;
        if (values[i] < property->minimum) property->minimum = values[i];
      }
    }
  }

  // Return minimum
  return minimum;
}



RNScalar R3MeshProperty::
Maximum(void) const
{
  // Compute maximum of property values
  if (maximum == RN_UNKNOWN) {
    R3MeshProperty *property = (R3MeshProperty *) this;
    property->maximum = 0;
    if (nvalues > 0) {
      property->maximum = -FLT_MAX;
      for (int i = 0; i < nvalues; i++) {
        if (values[i] == RN_UNKNOWN) continue;
        if (values[i] > property->maximum) property->maximum = values[i];
      }
    }
  }

  // Return maximum
  return maximum;
}



RNScalar R3MeshProperty::
StandardDeviation(void) const
{
  // Compute standard deviation of property values
  if (stddev == RN_UNKNOWN) {
    R3MeshProperty *property = (R3MeshProperty *) this;
    property->stddev = 0;
    if (nvalues > 0) {
      // Comute sum of squared residuals
      int count = 0;
      RNScalar ssd = 0;
      RNScalar avg = Mean();
      for (int i = 0; i < nvalues; i++) {
        if (values[i] == RN_UNKNOWN) continue;
        RNScalar delta = values[i] - avg;
        ssd += delta * delta;
        count++;
      }

      // Compute standard deviation
      if (count > 0) property->stddev = sqrt(ssd / count);
    }
  }

  // Return standard deviation
  return stddev;
}



RNScalar R3MeshProperty::
Variance(void) const
{
  // Return variance of property values
  return StandardDeviation() * StandardDeviation();
}



RNScalar R3MeshProperty::
Median(void) const
{
  // Compute mean of property values
  if (median == RN_UNKNOWN) {
    R3MeshProperty *property = (R3MeshProperty *) this;
    property->median = Percentile(50);
  }

  // Return median
  return median;
}



RNScalar R3MeshProperty::
Percentile(RNScalar percentile) const
{
  // Compute value at given percentile (100=max, 0=min)
  if (nvalues == 0) return 0;

  // Copy values
  int count = 0;
  RNScalar *copy = new RNScalar [ nvalues ];
  for (int i = 0; i < nvalues; i++) {
    if (values[i] == RN_UNKNOWN) continue;
    copy[count++] = values[i];
  }

  // Check count of values
  if (count == 0) return 0;

  // Compute value at given percentile
  qsort(copy, count, sizeof(RNScalar), RNCompareScalars);
  int index = (int) (count * percentile / 100.0);
  if (index >= count) index = count - 1;
  RNScalar result = copy[index];

  // Delete copy
  delete [] copy;

  // Return result
  return result;
}



RNScalar R3MeshProperty::
L2Norm(void) const
{
  // Compute L2 norm of property values
  if (l2norm == RN_UNKNOWN) {
    R3MeshProperty *property = (R3MeshProperty *) this;
    property->l2norm = 0;
    if (nvalues > 0) {
      // Compute sum 
      int count = 0;
      RNScalar sum = 0;
      for (int i = 0; i < nvalues; i++) {
        if (values[i] == RN_UNKNOWN) continue;
        sum += values[i] * values[i];
        count++;
      }

      // Compute l2norm
      if (count > 0) property->l2norm = sqrt( sum / count );
    }
  }

  // Return l2norm
  return l2norm;
}



////////////////////////////////////////////////////////////////////////
// Manipulation functions
////////////////////////////////////////////////////////////////////////

void R3MeshProperty::
Abs(void)
{
  // Replace every value by its absolute value
  for (int i = 0; i < nvalues; i++) {
    if (values[i] == RN_UNKNOWN) continue;
    values[i] = fabs(values[i]);
  }
}



void R3MeshProperty::
Sqrt(void)
{
  // Replace every value by its sqrt
  for (int i = 0; i < nvalues; i++) {
    if (values[i] == RN_UNKNOWN) continue;
    values[i] = sqrt(values[i]);
  }
}



void R3MeshProperty::
Square(void)
{
  // Replace every value by its square
  for (int i = 0; i < nvalues; i++) {
    if (values[i] == RN_UNKNOWN) continue;
    values[i] = values[i] * values[i];
  }
}



void R3MeshProperty::
Negate(void)
{
  // Negate every value
  for (int i = 0; i < nvalues; i++) {
    if (values[i] == RN_UNKNOWN) continue;
    values[i] = -values[i];
  }
}



void R3MeshProperty::
Invert(void)
{
  // Invert every value
  for (int i = 0; i < nvalues; i++) {
    if (values[i] == RN_UNKNOWN) continue;
    if (values[i] == 0) continue;
    values[i] = 1.0 / values[i];
  }
}



void R3MeshProperty::
Clear(RNScalar value)
{
  // Set every value
  for (int i = 0; i < nvalues; i++) {
    values[i] = value;
  }
}



void R3MeshProperty::
Copy(const R3MeshProperty& property)
{
  // Copy every value
  for (int i = 0; i < nvalues; i++) {
    values[i] = property.values[i];
  }
}



void R3MeshProperty::
Add(RNScalar value)
{
  // Add scalar to every value
  for (int i = 0; i < nvalues; i++) {
    if (values[i] == RN_UNKNOWN) continue;
    values[i] += value;
  }
}



void R3MeshProperty::
Add(const R3MeshProperty& property)
{
  // Add every value 
  for (int i = 0; i < nvalues; i++) {
    if (values[i] == RN_UNKNOWN) continue;
    if (property.values[i] == RN_UNKNOWN) continue;
    values[i] += property.values[i];
  }
}



void R3MeshProperty::
Subtract(RNScalar value)
{
  // Subtract scalar from every value
  for (int i = 0; i < nvalues; i++) {
    if (values[i] == RN_UNKNOWN) continue;
    values[i] -= value;
  }
}



void R3MeshProperty::
Subtract(const R3MeshProperty& property)
{
  // Subtract every value 
  for (int i = 0; i < nvalues; i++) {
    if (values[i] == RN_UNKNOWN) continue;
    if (property.values[i] == RN_UNKNOWN) continue;
    values[i] -= property.values[i];
  }
}



void R3MeshProperty::
Multiply(RNScalar value)
{
  // Multiply every value by a scalar
  for (int i = 0; i < nvalues; i++) {
    if (values[i] == RN_UNKNOWN) continue;
    values[i] *= value;
  }
}



void R3MeshProperty::
Multiply(const R3MeshProperty& property)
{
  // Multiply every value 
  for (int i = 0; i < nvalues; i++) {
    if (values[i] == RN_UNKNOWN) continue;
    if (property.values[i] == RN_UNKNOWN) continue;
    values[i] *= property.values[i];
  }
}



void R3MeshProperty::
Divide(RNScalar value)
{
  // Check value
  if (value == 0) return;

  // Divide every value by a scalar
  for (int i = 0; i < nvalues; i++) {
    if (values[i] == RN_UNKNOWN) continue;
    values[i] /= value;
  }
}



void R3MeshProperty::
Divide(const R3MeshProperty& property)
{
  // Divide every value 
  for (int i = 0; i < nvalues; i++) {
    if (values[i] == RN_UNKNOWN) continue;
    if (property.values[i] == RN_UNKNOWN) continue;
    if (property.values[i] == 0) continue;
    values[i] /= property.values[i];
  }
}



void R3MeshProperty::
Pow(RNScalar exponent)
{
  // Raise every value to an exponent
  for (int i = 0; i < nvalues; i++) {
    if (values[i] == RN_UNKNOWN) continue;
    values[i] = pow(values[i], exponent);
  }
}



void R3MeshProperty::
Threshold(RNScalar threshold, RNScalar low, RNScalar high)
{
  // Threshold every value
  for (int i = 0; i < nvalues; i++) {
    if (values[i] == RN_UNKNOWN) continue;
    if (values[i] <= threshold) values[i] = low;
    else values[i] = high;
  }
}



void R3MeshProperty::
AddVertexValue(int vertex_index, RNScalar value)
{
  // Add value to a vertex
  values[vertex_index] += value;
}



void R3MeshProperty::
AddVertexValue(R3MeshVertex *vertex, RNScalar value)
{
  // Add value to a vertex
  AddVertexValue(mesh->VertexID(vertex), value);
}



void R3MeshProperty::
SetVertexValue(int vertex_index, RNScalar value)
{
  // Set vertex value
  values[vertex_index] = value;
}



void R3MeshProperty::
SetVertexValue(R3MeshVertex *vertex, RNScalar value)
{
  // Set vertex value
  int vertex_index = mesh->VertexID(vertex);
  SetVertexValue(vertex_index, value);
}



void R3MeshProperty::
SetName(char *name)
{
  // Set name
  strncpy(this->name, name, 1024);
  this->name[0123] = '\0';
}



void R3MeshProperty::
Normalize(void)
{
  // Get mean and stddev
  RNScalar m = Mean();
  RNScalar s = StandardDeviation();
  if (s == 0) return;

  // Normalize by mean and stddev
  for (int i = 0; i < nvalues; i++) {
    if (values[i] == RN_UNKNOWN) continue;
    values[i] = (values[i] - m) / s;
  }
}



void R3MeshProperty::
Percentilize(void)
{
  // Make array of vertices
  R3MeshVertex **vertices = new R3MeshVertex *[ mesh->NVertices() ];
  for (int i = 0; i < mesh->NVertices(); i++) vertices[i] = mesh->Vertex(i);

  // Sort vertices array according to value
  for (int i1 = 0; i1 < mesh->NVertices(); i1++) {
    for (int i2 = i1+1; i2 < mesh->NVertices(); i2++) {
      if (values[mesh->VertexID(vertices[i2])] < values[mesh->VertexID(vertices[i1])]) {
        R3MeshVertex *swap = vertices[i1];
        vertices[i1] = vertices[i2];
        vertices[i2] = swap;
      }
    }
  }

  // Replace values with percentiles
  for (int i = 0; i < mesh->NVertices(); i++) {
    R3MeshVertex *vertex = vertices[i];
    int vertex_index = mesh->VertexID(vertex);
    values[vertex_index] = 100.0 * i / mesh->NVertices();
  }
}



void R3MeshProperty::
Blur(RNScalar sigma)
{
  // Blur across mesh

  // Get convenient variables
  if (sigma == 0) return;
  RNScalar radius = 3 * sigma;
  RNScalar denom = -2 * sigma * sigma;

  // Make copy of values
  RNScalar *old_values = new RNScalar [ mesh->NVertices() ];
  for (int i = 0; i < mesh->NVertices(); i++) old_values[i] = values[i];

  // Blur value at every vertex
  for (int i = 0; i < mesh->NVertices(); i++) {
    R3MeshVertex *vertex = mesh->Vertex(i);

    // Compute distances
    RNLength *distances = mesh->DijkstraDistances(vertex, radius);

    // Compute blurred value
    RNScalar total_value = 0;
    RNScalar total_weight = 0;
    for (int j = 0; j < mesh->NVertices(); j++) {
      R3MeshVertex *neighbor_vertex = mesh->Vertex(j);
      int neighbor_id = mesh->VertexID(neighbor_vertex);
      RNLength distance = distances[neighbor_id];
      if (distance > radius) continue;
      RNScalar value = old_values[neighbor_id];
      RNScalar weight = exp(distance * distance / denom); 
      total_value += weight * value;
      total_weight += weight;
    }

    // Assign blurred value (normalized by total weight)
    if (total_weight > 0) {
      values[i] = total_value / total_weight;
    }

    // Delete distances
    delete [] distances;
  }

  // Delete copy of values
  delete [] old_values;
}



void R3MeshProperty::
DoG(RNScalar sigma)
{
  // Difference of Gaussians
  // Replace values with difference from blurred version
  R3MeshProperty blur(*this); 
  blur.Blur(sigma);
  Subtract(blur);
}



void R3MeshProperty::
Strength(RNScalar sigma)
{
  // Replace values with significance score
  // For now, compute product of value and DoG
  R3MeshProperty dog(*this); 
  dog.DoG(sigma);
  dog.Abs();
  Multiply(dog);
}



// Include file for mesh property class



// Class definition

struct R3MeshProperty {
public:
  // Constructors/destructors
  R3MeshProperty(R3Mesh *mesh, const char *name = NULL, const RNScalar *vertex_values = NULL);
  R3MeshProperty(const R3MeshProperty& property);
  ~R3MeshProperty(void);

  // Property info functions
  R3Mesh *Mesh(void) const;
  const char *Name(void) const;

  // Value access functions
  int NVertexValues(void) const;
  RNScalar VertexValue(R3MeshVertex *vertex) const;
  RNScalar VertexValue(int vertex_index) const;

  // Statistics functions
  RNScalar Mean(void) const;
  RNScalar Variance(void) const;
  RNScalar StandardDeviation(void) const;
  RNScalar Minimum(void) const;
  RNScalar Maximum(void) const;
  RNScalar Median(void) const;
  RNScalar Percentile(RNScalar percentile) const;
  RNScalar L1Norm(void) const;
  RNScalar L2Norm(void) const;

  // Manipulation functions
  void Abs(void);
  void Sqrt(void);
  void Square(void);
  void Negate(void);
  void Invert(void);
  void Clear(RNScalar value = 0);
  void Copy(const R3MeshProperty& property);
  void Add(RNScalar value);
  void Add(const R3MeshProperty& property);
  void Subtract(RNScalar value);
  void Subtract(const R3MeshProperty& property);
  void Multiply(RNScalar value);
  void Multiply(const R3MeshProperty& property);
  void Divide(RNScalar value);
  void Divide(const R3MeshProperty& property);
  void Pow(RNScalar exponent);
  void Threshold(RNScalar threshold, RNScalar low, RNScalar high);
  void AddVertexValue(R3MeshVertex *vertex, RNScalar value);
  void AddVertexValue(int vertex_index, RNScalar value);
  void SetVertexValue(R3MeshVertex *vertex, RNScalar value);
  void SetVertexValue(int vertex_index, RNScalar value);
  void SetName(char *name);

  // More property manipulation functions
  void Normalize(void);
  void Percentilize(void);
  void Blur(RNScalar sigma);
  void DoG(RNScalar sigma);
  void Strength(RNScalar sigma);

  // Arithmetic operators
  R3MeshProperty& operator=(const R3MeshProperty& property);
  R3MeshProperty& operator+=(RNScalar scale);
  R3MeshProperty& operator+=(const R3MeshProperty& property);
  R3MeshProperty& operator-=(RNScalar scale);
  R3MeshProperty& operator-=(const R3MeshProperty& property);
  R3MeshProperty& operator*=(RNScalar scale);
  R3MeshProperty& operator*=(const R3MeshProperty& property);
  R3MeshProperty& operator/=(RNScalar scale);
  R3MeshProperty& operator/=(const R3MeshProperty& property);

public:
  R3Mesh *mesh;
  char name[1024];
  int nvalues;
  RNScalar *values;
  RNScalar mean;
  RNScalar stddev;
  RNScalar minimum;
  RNScalar maximum;
  RNScalar median;
  RNScalar l2norm;
};



// Inline functions

inline R3Mesh *R3MeshProperty::
Mesh(void) const
{
  // Return mesh
  return mesh;
}



inline const char *R3MeshProperty::
Name(void) const
{
  // Return property name
  return name;
}



inline int R3MeshProperty::
NVertexValues(void) const
{
  // Return number of vertex values
  return nvalues;
}




inline RNScalar R3MeshProperty::
VertexValue(int vertex_index) const
{
  // Return value at vertex
  return values[vertex_index];
}




inline RNScalar R3MeshProperty::
VertexValue(R3MeshVertex *vertex) const
{
  // Return value at vertex
  int vertex_index = mesh->VertexID(vertex);
  return VertexValue(vertex_index);
}




inline RNScalar R3MeshProperty::
L1Norm(void) const
{
  // Return L1Norm (mean) of property values
  return Mean();
}



inline R3MeshProperty& R3MeshProperty::
operator+=(RNScalar value) 
{
  // Add value to all property values 
  Add(value);
  return *this;
}



inline R3MeshProperty& R3MeshProperty::
operator+=(const R3MeshProperty& property) 
{
  // Add passed property values to corresponding entries of this property
  Add(property);
  return *this;
}



inline R3MeshProperty& R3MeshProperty::
operator-=(RNScalar value) 
{
  // Subtract value from all property values 
  Subtract(value);
  return *this;
}



inline R3MeshProperty& R3MeshProperty::
operator-=(const R3MeshProperty& property) 
{
  // Subtract passed property values from corresponding entries of this property
  Subtract(property);
  return *this;
}



inline R3MeshProperty& R3MeshProperty::
operator*=(RNScalar value) 
{
  // Multiply property values by value
  Multiply(value);
  return *this;
}



inline R3MeshProperty& R3MeshProperty::
operator*=(const R3MeshProperty& property) 
{
  // Multiply passed property values by corresponding entries of this property
  Multiply(property);
  return *this;
}



inline R3MeshProperty& R3MeshProperty::
operator/=(RNScalar value) 
{
  // Divide property values by value
  Divide(value);
  return *this;
}



inline R3MeshProperty& R3MeshProperty::
operator/=(const R3MeshProperty& property) 
{
  // Divide passed property values by corresponding entries of this property
  Divide(property);
  return *this;
}




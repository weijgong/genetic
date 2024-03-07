//------------------------------------------------------------------------------
//
// SAT_VecMat.h
// 
// Purpose: 
//
//   Vector/matrix operations
//   
// Last modified:
//
//   2000/03/04  OMO  Final version (1st edition)
//   2005/04/14  OMO  Final version (2nd reprint)
//
// (c) 1999-2005  O. Montenbruck, E. Gill
//
//------------------------------------------------------------------------------

#ifndef INC_SAT_VECMAT_H
#define INC_SAT_VECMAT_H

#include <iostream>


//------------------------------------------------------------------------------
//
// Vector (class definition)
//
// Purpose:
//
//   Vector data type and associated operations
//
//------------------------------------------------------------------------------
class Matrix;

class Vector
{

  public:
    
    friend class Matrix;
    
    // Constructors
    Vector ();                              // Vector without elements
    Vector (int Size);                      // Nullvector of specified size 
    Vector (const Vector& V);               // Vector copy
    Vector (const double* p, int N);        // Array copy
    Vector (double x, double y, double z);  // 3dim-Vector
    Vector (double x, double y, double z,   // 6dim-Vector
            double X, double Y, double Z);  
	Vector(double x, double y, double z,   // 7dim-Vector
		double X, double Y, double Z, double my);

    // Destructor
    ~Vector();
    
    // Size
    int size() const { return n; };
    Vector& resize(int Size);
    
    // Assignment
    Vector& operator=(const double value);
    Vector& operator=(const Vector& V);

    // Component access (Fortran notation)
    double  operator () (int i) const { return v[i]; };
    double& operator () (int i)       { return v[i]; };
    Vector slice (int first, int last) const;

    // Square root of vector elements
    Vector Sqrt() const;

    // Concatenation 
    friend Vector operator &(const Vector& a, double b);
    friend Vector operator &(double a, const Vector& b);
    friend Vector operator &(const Vector& a, const Vector& b);
    friend Vector Stack     (const Vector& a, const Vector& b);

    // Vector from polar angles
    friend Vector VecPolar (double azim, double elev, double r);

    // Vector addition/subtraction with assignment
    void operator += (const Vector& V);
    void operator -= (const Vector& V);
    
    // Dot product, norm, cross product
    friend double Dot (const Vector& left, const Vector& right);
	  friend double Ang (const Vector& left, const Vector& right);
    friend double Norm (const Vector& V);
    friend Vector Cross (const Vector& left, const Vector& right);

    // Scalar multiplication and division of a vector
    friend Vector operator * (double value, const Vector& V);
    friend Vector operator * (const Vector& V, double value);
    friend Vector operator / (const Vector& V, double value);

    // Negation of a vector (unary minus)
    friend Vector operator - (const Vector& V);

    // Vector addition and subtraction
    friend Vector operator + (const Vector& left, const Vector& right);    
    friend Vector operator - (const Vector& left, const Vector& right);    
    
    // Diagonal matrix
    friend Matrix Diag(const Vector& Vec);

    // Vector/matrix product
    friend Vector operator * (const Matrix& Mat, const Vector& Vec);
    friend Vector operator * (const Vector& Vec, const Matrix& Mat);

    // Dyadic product
    friend Matrix Dyadic (const Vector& left, const Vector& right);

    // Output
    friend std::ostream& operator << (std::ostream& os, const Vector& Vec);

  private:
        
    // Elements
    int    n;      // Dimension
    double *v;     // Vector v(n)

};


//------------------------------------------------------------------------------
//
// Matrix (class definition)
//
// Purpose:
//
//   Matrix data type and associated operations
//
//------------------------------------------------------------------------------

class Matrix
{

  public:

    // Constructors
    Matrix ();                                      // Matrix without elements
    Matrix (int dim1, int dim2);                    // Nullmatrix 
    Matrix (const Matrix& M_);                      // Matrix copy
    Matrix (const double* p, int dim1, int dim2);   // Array copy

    // Destructor
    ~Matrix();

    // Assignment
    Matrix& operator=(const double value);
    Matrix& operator=(const Matrix& M_);

    // Size
    int size1() const { return n; };
    int size2() const { return m; };
    Matrix& resize(int dim1, int dim2);
    
    // Component access (Fortran notation)
    double  operator () (int i, int j) const { return M[i][j]; };   
    double& operator () (int i, int j)       { return M[i][j]; };   
    Vector Col(int j) const;      
    Vector Row(int i) const;      
    Vector Diag() const;
    double Trace() const;
    double Trace(int low, int upp) const;
    Matrix slice(int first_row, int last_row, int first_col, int last_col);
    void SetCol(int j, const Vector& Col);
    void SetRow(int i, const Vector& Row);

    // Concatenation 
    friend Matrix operator &(const Matrix& A, const Vector& Row);
    friend Matrix operator &(const Vector& Row, const Matrix& A);
    friend Matrix operator &(const Matrix& A, const Matrix& B);
    friend Matrix operator |(const Matrix& A, const Vector& Col);
    friend Matrix operator |(const Vector& Col, const Matrix& A);
    friend Matrix operator |(const Matrix& A, const Matrix& B);

    // Matrix addition/subtraction with assignment
    void operator += (const Matrix& V);
    void operator -= (const Matrix& V);

    // Unit matrix
    friend Matrix Id(int Size);

    // Diagonal matrix
    friend Matrix Diag(const Vector& Vec);

    // Elementary rotations
    friend Matrix R_x(double Angle);
    friend Matrix R_y(double Angle);
    friend Matrix R_z(double Angle);
	friend Matrix ECIOrbitMatrix(double ui, double in, double om);

    // Transposition and inverse
    friend Matrix Transp(const Matrix& Mat);
    friend Matrix Inv(const Matrix& Mat);

    // Scalar multiplication and division of a matrix
    friend Matrix operator * (double value, const Matrix& Mat);
    friend Matrix operator * (const Matrix& Mat, double value);
    friend Matrix operator / (const Matrix& Mat, double value);

    // Unary minus
    friend Matrix operator - (const Matrix& Mat);

    // Matrix addition and subtraction
    friend Matrix operator + (const Matrix& left, const Matrix& right);
    friend Matrix operator - (const Matrix& left, const Matrix& right);    

    // Matrix product
    friend Matrix operator * (const Matrix& left, const Matrix& right);
    
    // Vector/matrix product
    friend Vector operator * (const Matrix& Mat, const Vector& Vec);
    friend Vector operator * (const Vector& Vec, const Matrix& Mat);

    // Dyadic product
    friend Matrix Dyadic (const Vector& left, const Vector& right);

    // Output
    friend std::ostream& operator << (std::ostream& os, const Matrix& Mat);

  private:

    // Elements
    int      n;                       // First dimension (number of rows)
    int      m;                       // Second dimension (number of columns)
    double **M;                       // Matrix M(n,m)

};


//------------------------------------------------------------------------------
//
// Basic Linear Algebra
//
//------------------------------------------------------------------------------

void LU_Decomp ( Matrix& A, Vector& Indx );
void LU_BackSub ( Matrix& A, Vector& Indx, Vector& b );


#endif  // include-Blocker

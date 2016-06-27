/**
* @file    Vertex.hpp
* @author  Marcos Lage      <mlage@mat.puc-rio.br>
* @author  Math Dept, PUC-Rio
* @author  Lab Matm�dia
* @date    14/02/2006
*
* @brief  (Vertex Class).
*/
/**------------------------------------------------------------------------------------*/
#ifndef _VERTEX_HPP_
#define _VERTEX_HPP_

#include "Defs.hpp"

/** Vertex id type; */
typedef int Vid;

using namespace std;
/**------------------------------------------------------------------------------------*/
/** Vertex class for CHE data structure 
  * \brief Vertex class*/
class Vertex
/**------------------------------------------------------------------------------------*/
{
//-- Vertex protected data --//
protected:
  /** Vertex coordenates
    * \brief coordenate x of the vertex*/
  float  _x;    
  /** Vertex coordenates
    * \brief coordenate y of the vertex*/
  float  _y;  
  /** Vertex coordenates
    * \brief coordenate z of the vertex*/
  float  _z;
  /** Vertex normal coordenates
    * \brief coordenate x of the vertex normal*/
  float _nx;
   /** Vertex normal coordenates
    * \brief coordenate y of the vertex normal*/
  float _ny;
  /** Vertex normal coordenates
    * \brief coordenate z of the vertex normal*/
  float _nz;
  
/*--------------------------------------------------------------------------armaz lap das normais -------------------------------------------*/
	float _lnx;
	/** Vertex normal coordenates (Laplaciano)				
	 * \brief coordenate y of the vertex normal*/
	float _lny;
	/** Vertex normal coordenates(Laplaciano)
	 * \brief coordenate z of the vertex normal*/
	float _lnz;
/*--------------------------------------------------------------------------------------------------------------------------------------------*/
	
  	float _grad1;
	float _grad2;
	float _grad3;
	
	
  /** Property Dimesion
  * \brief dimension of vector of properties*/
  int    _nf;
  /** Vertex vector of properties
    * \brief Property of each vertex.*/
  float* _f;
	

	


/**------------------------------------------------------------------------------------*/
public:
/**------------------------------------------------------------------------------------*/
  /** \brief Default constructor.*/
	Vertex():_x(0), _y(0), _z(0), _nx(0), _ny(0), _nz(0),  _lnx(0), _lny(0), _lnz(0),_grad1(0), _grad2(0), _grad3(0),_nf(0),_f(NULL) {} 

  /** \brief First constructor*/
	Vertex( float  x, float  y, float  z ):_x(x), _y(y), _z(z), _nx(0), _ny(0), _nz(0),_lnx(0), _lny(0),_lnz(0),_grad1(0), _grad2(0), _grad3(0), _nf(0) , _f(NULL){ } 

  /** \brief Second constructor*/  
	Vertex( float  x, float  y, float  z,
		    float nx, float ny, float nz,
			float lnx, float lny, float lnz,
		    float grad1, float grad2, float grad3,
			float* f, int nf ): _x(x), _y(y), _z(z), _nx(nx), _ny(ny), _nz(nz),_lnx(lnx), _lny(lny), _lnz(lnz),_grad1(grad1), _grad2(grad2), _grad3(grad3), _f(f),_nf(nf) {} 

   /** \brief Copy constructor
     * \param Vertex& v.*/
	Vertex(const Vertex& v):_x(v.x()), _y(v.y()), _z(v.z()), _nx(v.nx()), _ny(v.ny()), _nz(v.nz()),_lnx(v.lnx()), _lny(v.lny()), _lnz(v.lnz()),_grad1(v.grad1()), _grad2(v.grad2()), _grad3(v.grad3()), _f(v.f()),_nf(v.nf()) {}	

   /** \brief Destructor.*/
	~Vertex(){ erase_f(); }

   /** \brief Assignment operator.*/
    Vertex& operator = (const Vertex& v){
		_x = v.x(); _y = v.y(); _z = v.z(); _nx = v.nx(); _ny = v.ny(); _nz = v.nz();_lnx = v.lnx(); _lny = v.lny(); _lnz = v.lnz(); _grad1=v.grad1(); _grad2=v.grad2(); _grad3=v.grad3();   _nf = v._nf; _f = v.f();		
	  return *this;
   }
/**------------------------------------------------------------------------------------*/
//-- Vertex methods. --//
/**------------------------------------------------------------------------------------*/
public:
	/** \brief Access to the coordinate x of the vertex*/
	inline const float     x() const{ return _x ; }

	/** \brief Access to the coordinate y of the vertex*/
	inline const float     y() const{ return _y ; }

	/** \brief Access to the coordinate z of the vertex*/
	inline const float     z() const{ return _z ; }

	/** \brief Access to the coordenate nx of the normal*/
	inline const float    nx() const{ return _nx; } 

	/** \brief Access to the coordenate ny of the normal*/
	inline const float    ny() const{ return _ny; }

	/** \brief Access to the coordenate nz of the normal*/
	inline const float    nz() const{ return _nz; }
	
	
	
	/** \brief Access to the coordenate lnx of the normal (laplaciano)*/		
	inline const float    lnx() const{ return _lnx; } 
	
	/** \brief Access to the coordenate lny of the normal(laplaciano)*/
	inline const float    lny() const{ return _lny; }
	
	/** \brief Access to the coordenate lnz of the normal(laplaciano)*/
	inline const float    lnz() const{ return _lnz; }
	

	inline const float    grad1() const{ return _grad1; } 
	inline const float    grad2() const{ return _grad2; } 
	inline const float    grad3() const{ return _grad3; } 
	
	
	
	
	
	/** \brief Access to the number of properties*/
	inline const int  nf() const{ return _nf; }

	/** \brief Access to the coordenate f of the field*/
	inline float* f() const{ return  _f; } 
/**------------------------------------------------------------------------------------*/
public:
	/** \brief Sets the coordenate x of the vertex
	  * \param const float x */
	inline const void  set_x (const float x ) { _x=x  ; }

	/** \brief Sets the coordenate y of the vertex
	* \param const float y */
	inline const void  set_y (const float y ) { _y=y  ; }

	/** \brief Sets the coordenate z of the vertex
	  * \param const float z */
	inline const void  set_z (const float z ) { _z=z  ; }

	/** \brief Sets the coordenate nx of the normal
	  * \param const float nx */
	inline const void  set_nx(const float nx) { _nx=nx; }

	/** \brief Sets the coordenate ny of the normal
	* \param const float ny */
	inline const void  set_ny(const float ny) { _ny=ny; }

	/** \brief Sets the coordenate nz of the normal
	  * \param const float nz */
	inline const void  set_nz(const float nz) { _nz=nz; }
	
	
	/** \brief Sets the coordenate lnx of the normal (laplaciano)		
	 * \param const float lnx */
	inline  const void  set_lnx(const float lnx) { _lnx=lnx; }
	
	/** \brief Sets the coordenate lny of the normal (laplaciano)
	 * \param const float lny */
	inline const void  set_lny(const float lny) { _lny=lny; }
	
	/** \brief Sets the coordenate lnz of the normal (laplaciano)
	 * \param const float lnz */
	inline const void  set_lnz(const float lnz) { _lnz=lnz; }
	
	inline const void  set_grad1(const float grad1) { _grad1=grad1; }
	
	inline const void  set_grad2(const float grad2) { _grad2=grad2; }
	
	inline const void  set_grad3(const float grad3) { _grad3=grad3; }

		
	/** \brief acrescenta um valor a coordenada x do laplaciano da normal		
	 * \param const float lnx */
	inline  const void  add_lnx(const float lnx) { _lnx+=lnx; }
	
	/** \brief acrescenta um valor a coordenada y do laplaciano da normal
	 * \param const float lny */
	inline const void  add_lny(const float lny) { _lny+=lny; }
	
	/** \brief acrescenta um valor a coordenada z do laplaciano da normal
	 * \param const float lnz */
	inline const void  add_lnz(const float lnz) { _lnz+=lnz; }
	
	inline const void  add_grad1(const float grad1) { _grad1+=grad1; }
	
	inline const void  add_grad2(const float grad2) { _grad2+=grad2; }
	
	inline const void  add_grad3(const float grad3) { _grad3+=grad3; }
	
	
	/** \brief multiplica um valor a coordenada x do laplaciano da normal		
	 * \param const float lnx */
	inline  const void  mult_lnx(const float lnx) { _lnx*=lnx; }
	
	/** \brief multiplica um valor a coordenada y do laplaciano da normal
	 * \param const float lny */
	inline const void  mult_lny(const float lny) { _lny*=lny; }
	
	/** \brief multiplica um valor a coordenada z do laplaciano da normal
	 * \param const float lnz */
	inline const void  mult_lnz(const float lnz) { _lnz*=lnz; }
	
	

	
	
	/** \brief Sets the coordenate nz of the normal
	 * \param const float nz */
	inline const void  set_nf(const int nf) { _nf=nf; }
	
	/** \brief Sets the coordenate fx of the field
	 * \param const float fx */
	inline const void  set_f(float* f) { _f=f; }
	
	/** \brief Sets the coordenate fx of the field
	 * \param const float fx */
	inline const void erase_f() { if(_f && _nf == 0) { delete _f; } if(_f && _nf > 0) { delete [] _f; } }
	
	/** \brief Sets the coordenate fx of the field
	 * \param const float fx */
	inline const void alloc_f(int _size) { erase_f(); _f = new float[_size]; } 

/**------------------------------------------------------------------------------------*/
//-- Geometric Operations --//
/**------------------------------------------------------------------------------------*/
public:
	/** \brief Computes the normal of a triangle
    * \param const Vertex &v0
    * \param const Vertex &v1
    * \param const Vertex &v2
    * \param float* n*/
  static void normal ( const Vertex &v0, const Vertex &v1, const Vertex &v2, float *n)
  {
      float a0[3], a1[3];

      a0[0]= v1.x()-v0.x();
      a0[1]= v1.y()-v0.y();
      a0[2]= v1.z()-v0.z();

      a1[0]= v2.x()-v0.x();
      a1[1]= v2.y()-v0.y();
      a1[2]= v2.z()-v0.z();

      n[0]= a0[1]*a1[2]-a1[1]*a0[2];
      n[1]= a0[2]*a1[0]-a1[2]*a0[0];
      n[2]= a0[0]*a1[1]-a1[0]*a0[1];
	  
	  normalize(n);
  }
/**------------------------------------------------------------------------------------*/
  /** \brief computes the norm of a vector in R^3
    * \param const float a
    * \param const float b
    * \param const float c*/
  inline static float norm  ( const float a, const float b, const float c ) { return sqrt( a*a + b*b + c*c ) ; }

  /** \brief Normalizes a vector in R^3
    * \param float* v*/
  static void normalize( float *v )
  {
    float n = norm (v[0], v[1], v[2]);
    if(n > FLT_EPSILON )
    {
      v[0]/=n;
      v[1]/=n;
      v[2]/=n;
    }
  }
/**------------------------------------------------------------------------------------*/
 
    /** \brief Vertex minus operation */
	friend Vertex operator - ( Vertex& vf,  Vertex& v0)
	{
		Vertex p(vf);

		p._x -= v0._x;
		p._y -= v0._y;
		p._z -= v0._z;

		p._nx -= v0._nx;
		p._ny -= v0._ny;
		p._nz -= v0._nz;

		if(vf._nf != v0._nf) { cout << "Vertex: Minus op ERROR" << endl; exit(1);}

		for(int i=0; i<vf._nf; ++i)
			p._f[i] -= v0._f[i];

		return p;
	}

/**------------------------------------------------------------------------------------*/

    /** \brief Vertex minus operation */
	friend Vertex operator + ( Vertex& vf,  Vertex& v0)
	{
		Vertex p(vf);

		p._x += v0._x;
		p._y += v0._y;
		p._z += v0._z;

		p._nx += v0._nx;
		p._ny += v0._ny;
		p._nz += v0._nz;

		if(vf._nf != v0._nf) { cout << "Vertex: Plus op ERROR" << endl; exit(1);}

		for(int i=0; i<vf._nf; ++i)
			p._f[i] += v0._f[i];

		return p;
	}
/**------------------------------------------------------------------------------------*/

	/** \brief Vertex Dot product */
	friend float operator * ( Vertex& v0,  Vertex& v1)
	{
		return (v0._x * v1._x + v0._y * v1._y + v0._z * v1._z);
	}

/**------------------------------------------------------------------------------------*/

	/** \brief Vertex Dot product */
	friend Vertex operator * ( float c, Vertex& v0 )
	{
		Vertex p(v0);

		p._x *= c;
		p._y *= c;
		p._z *= c;

		p._nx *= c;
		p._ny *= c;
		p._nz *= c;

		for(int i=0; i<v0._nf; ++i)
			p._f[i] *= c;

		return p;
	}

/**------------------------------------------------------------------------------------*/
  /** \brief Vertex cross product */
	friend Vertex operator ^ ( Vertex& v0,  Vertex& v1)
	{
		Vertex p;
		p._x = v0._y * v1._z - v0._z * v1._y ;
		p._y = v0._z * v1._x - v0._x * v1._z ;
		p._z = v0._x * v1._y - v0._y * v1._x ;
		return p;
	}

/**------------------------------------------------------------------------------------*/


/**------------------------------------------------------------------------------------*/
//-- Vertex bool operators --//
/**------------------------------------------------------------------------------------*/
public: 
  /** \brief Equal operator for Vertex objects 
    * \param Vertex& v */
   bool operator == (const Vertex& v){ return (x()==v.x() && y()==v.y() && z()==v.z() && nx()==v.nx() && ny()==v.ny() && nz()==v.nz() && f()==v.f() && nf()==v.nf()); }
  
  /** \brief Equal operator for Vertex objects 
    * \param Vertex& v */
   bool operator != (const Vertex& v){ return (x()!=v.x() || y()!=v.y() || z()!=v.z() || nx()!=v.nx() || ny()!=v.ny() || nz()!=v.nz() || f()!=v.f() || nf()!=v.nf());}
};
#endif
/**------------------------------------------------------------------------------------*/

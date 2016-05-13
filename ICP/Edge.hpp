/**
* @file    Edge.hpp
* @author  Marcos Lage      <mlage@mat.puc-rio.br>
* @author  Hélio  Lopes     <lopes@mat.puc-rio.br>
* @author  Math Dept, PUC-Rio
* @author  Lab Matmídia
* @version 
* @date     06/06/2005
*
* @brief  (Edge Class).
*/
/**------------------------------------------------------------------------------------*/

#ifndef _Edge_HPP_
#define _Edge_HPP_

#include "Defs.hpp"

using namespace std;
/**------------------------------------------------------------------------------------*/
/** Edge class for CHF data structure 
  * \brief Edge class*/
class Edge
/**------------------------------------------------------------------------------------*/
{
//-- Edge protected data --//
protected:
  /** Edge Vertices
    * \brief first vertice of the edge*/
	int  _a;    
  /** Edge Vertices
    * \brief second vertice of the edge*/
	int  _b;    
  /** Edge Property
    * \brief User defined property*/
	float _p;    

/**------------------------------------------------------------------------------------*/
//-- Edge constructors.--// 
/**------------------------------------------------------------------------------------*/
public:
  /** \brief Default constructor.*/
	Edge():_a(-1), _b(-1), _p(FLT_MAX){}

  /** \brief First constructor
    * \param int  a  Edge _a
    * \param int  b  Edge _b
    * \param int  h  Edge _h */
	Edge( int  a, int  b ):_a(a), _b(b), _p(FLT_MAX){}
  
  /** \brief Copy constructor
    * \param Edge& _e*/
	Edge( const Edge& e ):_a(e.a()), _b(e.b()), _p(e.p()){}

  /** \brief Destructor.*/
    ~Edge(){}

  /** \brief Assignment operator.*/
    Edge& operator = (const Edge& e){
		_a = e.a(); _b = e.b(); _p = e.p();    
		return *this;
	}
/**------------------------------------------------------------------------------------*/
//-- Edge methods.--//
/**------------------------------------------------------------------------------------*/
public:
	/** \brief Access to the vertex a of the edge*/
  inline const int   a() const{ return _a ; }

	/** \brief Access to the vertex b of the edge*/
  inline const int   b() const{ return _b ; }

	/** \brief Access to the a half-edge of the edge*/
  inline const float p() const{ return _p ; }

/**------------------------------------------------------------------------------------*/
public:
/**------------------------------------------------------------------------------------*/
  /** \brief Sets the vertex a of the edge
    * \param const int a */
  inline const void  set_a (const int a) { _a=a; }
	
  /** \brief Sets the vertex b of the edge
    * \param const int b */
  inline const void  set_b (const int b) { _b=b; }
	
  /** \brief Sets the a half-edge of the edge
    * \param const int h */
  inline const void  set_p (const float p) { _p=p  ; }

/**------------------------------------------------------------------------------------*/
//-- Edge i/o operators --//
/**------------------------------------------------------------------------------------*/
public:
  /** \brief Read operator for Edge objects 
    * \param istream& s
    * \param Edge& e */
  friend istream& operator >> (istream& s,  Edge& e)
  {
    int a, b; float p;

    s >> a >> b >> p  ;

    e.set_a(a); e.set_b(b); e.set_p(p);
	  
	return s;
  }
/**------------------------------------------------------------------------------------*/

  /** \brief Write operator for Edge objects 
    * \param ostream& s
    * \param Edge& e */
  friend ostream& operator << (ostream& s, const Edge& e)
  {
    s << "  Verts: ( " << e.a() << " , " << e.b() <<  " ) " << endl ;
    s << "   Prop:   " << e.p() << endl ;
	return s;
  }
/**------------------------------------------------------------------------------------*/

/**------------------------------------------------------------------------------------*/
//-- Edge bool operators --//
/**------------------------------------------------------------------------------------*/
public: 
  /** \brief Equal operator for Edge objects 
    * \param Edge& e */
   bool operator == (const Edge& e){ return (a()==e.a() && b()==e.b() && p()==e.p()) ;}
  
  /** \brief Equal operator for Edge objects 
    * \param Edge& e */
   bool operator != (const Edge& e){ return (a()!=e.a() && b()!=e.b() && p()!=e.p()) ;}
};
#endif
/**------------------------------------------------------------------------------------*/

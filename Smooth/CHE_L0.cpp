/**
* @file    CHE_L0.cpp
* @author  Marcos Lage      <mlage@mat.puc-rio.br>
* @author  Thomas Lewiner   <thomas.lewiner@polytechnique.org>
* @author  H�lio  Lopes     <lopes@mat.puc-rio.br>
* @author  Math Dept, PUC-Rio
* @author  Lab Matm�dia
* @date    14/02/2006
*
* @brief  (Compact Half-Edge Structure - Level 0)
*/
/*--------------------------------------------------------------------------*/
#include  "CHE_L0.hpp"
#include <GL/glut.h>
 
//#include <gl/glut.h>

using namespace std;
/*--------------------------------------------------------------------------*/
vector<Vid> CHE_L0::R_00(const Vid v)
//--------------------------------------------------//
/** Computes the vertices in the star of a given vertex.*/
{
  set<Vid>   sstar;
  vector<Vid> star;

  if( !v_valid(v) ) { 
    if(DBG) cout << "CHE_L0::Vertex Star ERROR: invalid vertex id" << endl;
    star.push_back(INV); 
    return star; 
  }

  for(HEid i=0; i<3*ntrig(); ++i)
    if( V(i) == v ) {
      sstar.insert( V(next(i)) );
      sstar.insert( V(prev(i)) );
    }

  // Copies in star
  for(set<Vid>::const_iterator iter = sstar.begin(); iter != sstar.end(); ++iter )
	star.push_back( *iter );

  sstar.clear();
  return star;
}
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
vector<TRid> CHE_L0::R_02(const Vid v)
//--------------------------------------------------//
/** Computes the triangles in the star of a given vertex.*/
{
  vector<Vid> star;

  if( !v_valid(v) ) 
  { 
    if(DBG) cout << "CHE_L0::Vertex Star ERROR: invalid vertex id" << endl;
    star.push_back(INV); 
    return star; 
  }

  for(HEid i=0; i<3*ntrig(); ++i)
    if( V(i) == v ) 
      star.push_back( trig(i) );

  return star;
}
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
vector<Vid> CHE_L0::R_10(const HEid h)
//--------------------------------------------------//
/** Computes the vertices in the star of a given edge.*/
{
  vector<Vid> star;

  if( !he_valid(h) ) 
  { 
    if(DBG) cout << "CHE_L0::Edge Star ERROR: invalid edge id" << endl;
    star.push_back(INV); 
    return star; 
  }

  Vid a = V(h); Vid b = V(next(h));
  star.push_back( V(prev(h)) );

  for(HEid i=0; i<3*ntrig(); ++i)
  { 
    if( trig(i) == trig(h) ) continue;

    if( (b == V(i) && a == V(next(i))) || (a == V(i) && b == V(next(i))) ){
      star.push_back( V(prev(i)) );
      break;
    }
  }
  return star;
}
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
vector<TRid> CHE_L0::R_12(const HEid h)
//--------------------------------------------------//
/** Computes the triangles in the star of a given edge.*/
{
  vector<Vid> star;

  if( !he_valid(h) ) 
  { 
    if(DBG) cout << "CHE_L0::Edge Star ERROR: invalid edge id" << endl;
    star.push_back(INV); 
    return star; 
  }

  // gets the edge's vertices 
  Vid a = V(h); Vid b = V(next(h));
  
  // push_back the first triangle
  star.push_back( trig(h) );

  for(HEid i=0; i<3*ntrig(); ++i)
  {
    // avoids triangle repetition
	if( trig(i) == trig(h) ) continue;
	
	// find the other half-edge (possible non-oriented)
    if( (b == V(i) && a == V(next(i))) ||  
		(a == V(i) && b == V(next(i))) ){
      star.push_back( trig(i) );
      break;
    }
  }
  return star;
}
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
vector<TRid> CHE_L0::R_22(const TRid t)
//--------------------------------------------------//
/** Computes the triangle incidents of a given triangle.*/
{
  vector<TRid> star;

  if( !tr_valid(t) ) 
  { 
    if(DBG) cout << "CHE_L0::Triangle Star ERROR: invalid triangle id" << endl;
    star.push_back(INV); 
    return star; 
  }

  Vid   a = V( 3*t ), 
        b = V(3*t+1),
        c = V(3*t+2);

  for(HEid i=0; i<3*ntrig(); ++i)
  {
	if( trig(i) == t ) continue;

    if( ( b == V(i) && a == V(next(i)) ) ||   
        ( a == V(i) && c == V(next(i)) ) ||   
        ( c == V(i) && b == V(next(i)) ) 
		||
		( a == V(i) && b == V(next(i)) ) ||
		( b == V(i) && c == V(next(i)) ) ||
		( c == V(i) && a == V(next(i)) ) )

      star.push_back( trig(i) );
  }
  return star;
}
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
void CHE_L0::compute_normals()
//--------------------------------------------------//
/** Computes the normals of the triangles.*/				//<---------------- nao eh a normal dos triangulos e sim dos vertices, ne?
{
  float norm[3];
  float  nrm[3];

  cout << "CHE_L0::compute_normals..." ;
/*
  for(Vid i=0; i< ntrig(); i++)
  {
    if( !v_valid(V(3*i)) || !v_valid(V(3*i+1)) || !v_valid(V(3*i+2)) ) continue;

    Vertex &v0 = G( V( 3*i ) );
    Vertex &v1 = G( V(3*i+1) );
    Vertex &v2 = G( V(3*i+2) );
*/  
	  for(int r=0; r< 3*ntrig(); r+=3)
	  {
		  
		  if( !tr_valid( trig(r) ) ) continue;
		  
		  Vertex &v0 = G(V(r));
		  Vertex &v1 = G(V(r+1));
		  Vertex &v2 = G(V(r+2));
		

    Vertex::normal( v0,v1,v2, norm ) ;
    Vertex::normalize(norm) ;
    //--- V0 ----//
	  nrm[0] = norm[0] + v0.nx();
	  nrm[1] = norm[1] + v0.ny();
	  nrm[2] = norm[2] + v0.nz();

    v0.set_nx( nrm[0] );
    v0.set_ny( nrm[1] );
    v0.set_nz( nrm[2] );					//acho que tem que normalizar o v0 aqui

    //--- V1 ----//
	  
//	Vertex::normalize(norm) ;
    nrm[0] = norm[0] + v1.nx();
    nrm[1] = norm[1] + v1.ny();
    nrm[2] = norm[2] + v1.nz();

            
    v1.set_nx( nrm[0] );
    v1.set_ny( nrm[1] );
    v1.set_nz( nrm[2] );

    //--- V2 ----//
//	Vertex::normalize(norm) ;
	  nrm[0] = norm[0] + v2.nx();
	  nrm[1] = norm[1] + v2.ny();
	  nrm[2] = norm[2] + v2.nz();

    

    v2.set_nx( nrm[0] );
    v2.set_ny( nrm[1] );
    v2.set_nz( nrm[2] );
  }
	flip_normals();
	normaliza_normais();
	
  cout << " done." << endl;
}
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
void CHE_L0::normaliza_normais()
/*--------------------------------------------------------------------------*/
{
	for (int i = 0 ; i < _nvert; i++){
		Vertex &P= G(i);  
		float N[3];
		
		N[0]=P.nx();
		N[1]=P.ny();
		N[2]=P.nz();
		
		Vertex::normalize(N);
		
		P.set_nx( N[0]);
		P.set_ny( N[1]);
		P.set_nz( N[2]);
	}
}
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
void CHE_L0::bounding_box( float *min, float *max )
//--------------------------------------------------//
{
  float t_mx,t_Mx,t_my,t_My,t_mz,t_Mz;

  t_mx=t_Mx=_G[0].x();
  t_my=t_My=_G[0].y();
  t_mz=t_Mz=_G[0].z();

  for(Vid i=1; i<nvert(); ++i)
  {
    if(_G[i].x() < t_mx ) t_mx=_G[i].x();
    if(_G[i].x() > t_Mx ) t_Mx=_G[i].x();
																
    if(_G[i].y() < t_my ) t_my=_G[i].y();
    if(_G[i].y() > t_My ) t_My=_G[i].y();
																
    if(_G[i].z() < t_mz ) t_mz=_G[i].z();
    if(_G[i].z() > t_Mz ) t_Mz=_G[i].z();
  }

  min[0]=t_mx; min[1]=t_my; min[2]=t_mz;
  max[0]=t_Mx; max[1]=t_My; max[2]=t_Mz;
}
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
void CHE_L0::legalize_model (float *min, float *max)
//--------------------------------------------------//
{
  float c[3];
  float l[3];
  float size=0;

  c[0]=(float)(max[0]+min[0])/2;
  c[1]=(float)(max[1]+min[1])/2;
  c[2]=(float)(max[2]+min[2])/2;

  for(Vid i=0; i<3; i++) {
    l[i]=fabs(max[i]-min[i]);
    if(l[i]>size) size=l[i];
  }

  for(Vid i=0; i<nvert(); i++)
  {
    float tx = _G[i].x();
    float ty = _G[i].y();
    float tz = _G[i].z();

     tx -= c[0];
    if(size != 0) tx /= 0.05*size;
    _G[i].set_x(tx);

     ty -= c[1];
    if(size != 0) ty /= 0.05*size;
    _G[i].set_y(ty);

     tz -= c[2];
    if(size != 0) tz /= 0.05*size;
    _G[i].set_z(tz);
  }

  for(Vid i=0; i<3; i++)
  {
    min[i] = -1.;
    max[i] =  1.;
  }
}
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
void CHE_L0::check()
//--------------------------------------------------//
/** Checks the mesh.*/
{
  if(   nvert() != static_cast<int>( _G.size()) ){ cout << "CHE_L0:: Erro nvert()!= G.size"   << endl; return;}
  if( 3*ntrig() != static_cast<int>( _V.size()) ){ cout << "CHE_L0:: Erro 3*ntrig()!= V.size" << endl; return;}

  for(int i=0; i<3*ntrig(); ++i)
  {
	if( !he_valid(i) ) continue;

    if( V(i) >= nvert() )
    {
      cout << "CHE_L0:: Erro V(" << i <<") >= nvert." << endl;
      return;
    }
    if( V(i) < 0 )
    {
      cout << "CHE_L0:: Erro V(" << i <<") < 0." << endl;
      return;
    }
  }
}
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
void CHE_L0::draw_smooth()
//--------------------------------------------------//
/** Draws the smooth surface with opengl.*/
{
	for(int r=0; r< 3*ntrig(); r+=3)
	{
		if( !tr_valid( trig(r) ) ) { continue; }

		const Vertex &v1 = G( V(r) );
		const Vertex &v2 = G(V(r+1));
		const Vertex &v3 = G(V(r+2));

		//glColor3f(0.2567,0.5,0.9);

		glBegin(GL_TRIANGLES);
			glNormal3d(v1.nx(),  v1.ny(),  v1.nz());
			glVertex3d(v1.x()  , v1.y(),   v1.z() );
							
			glNormal3d(v2.nx(),  v2.ny(),  v2.nz());
			glVertex3d(v2.x()  , v2.y(),   v2.z() );
						
			glNormal3d(v3.nx(),  v3.ny(),  v3.nz());
			glVertex3d(v3.x()  , v3.y(),   v3.z() );
		glEnd();
	}	
}
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
void CHE_L0::draw_wire()
//--------------------------------------------------//
/** Draws the surface in wireframe with opengl.*/
{
  for(int r=0; r< 3*ntrig(); r+=3)
  {
	if( !tr_valid( trig(r) ) ) continue;

	const Vertex &v1 = G(V(r));
	const Vertex &v2 = G(V(r+1));
	const Vertex &v3 = G(V(r+2));

	//glColor3f(0.65,0.65,0.65);
    glBegin( GL_LINES );
	    //---- e0 ----//
      glNormal3d(v1.nx(),  v1.ny(),  v1.nz());
      glVertex3d(v1.x() ,  v1.y() ,  v1.z() );
        
	    glNormal3d(v2.nx(),  v2.ny(),  v2.nz());
	    glVertex3d(v2.x() ,  v2.y() ,  v2.z() );
	    //---- e1 ----//
	    glNormal3d(v2.nx(),  v2.ny(),  v2.nz());
 	    glVertex3d(v2.x() ,  v2.y() ,  v2.z() );
                                
      glNormal3d(v3.nx(),  v3.ny(),  v3.nz());
	    glVertex3d(v3.x() ,  v3.y() ,  v3.z() );
	    //---- e2 ----//
	    glNormal3d(v3.nx(),  v3.ny(),  v3.nz());
	    glVertex3d(v3.x() ,  v3.y() ,  v3.z() );
                                         
      glNormal3d(v1.nx(),  v1.ny(),  v1.nz());
	    glVertex3d(v1.x() ,  v1.y() ,  v1.z() );
    glEnd();
  }
}
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
void CHE_L0::draw_verts()
//--------------------------------------------------//
/** Draws the verices of the surface with opengl.*/
{
	for(int i=0; i<nvert(); ++i)
	{
		if( !v_valid( i ) ) continue;
		
		const Vertex &v1 = G(i);
		glColor3f(0.65,0.65,0.65);
		glPointSize(3.5);
		glBegin(GL_POINTS);
			glNormal3d(v1.nx(),  v1.ny(),  v1.nz());
			glVertex3d(v1.x() ,  v1.y() ,  v1.z() );
		glEnd();
		glPointSize(1.0);
	}
}
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
void CHE_L0::flip_normals()
/*--------------------------------------------------------------------------*/
{
	for (int i = 0 ; i < _nvert; i++){
		float nx = _G[i].nx();
		float ny = _G[i].ny();
		float nz = _G[i].nz();
		
		_G[i].set_nx( -nx );
		_G[i].set_ny( -ny );
		_G[i].set_nz( -nz );
	}
}
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
void CHE_L0::draw_normals()
/*--------------------------------------------------------------------------*/
{
	float size = 1;
	float vert[3];
	float nrm[3];

	glColor3d(0.5,0.5,0.5);		
	
	for (int i = 0 ; i < _nvert; i++){
		vert[0] = _G[i].x();
		vert[1] = _G[i].y();
		vert[2] = _G[i].z();

		 nrm[0] = _G[i].nx();
		 nrm[1] = _G[i].ny();
		 nrm[2] = _G[i].nz();

		float dv = sqrt(nrm[0]*nrm[0] + nrm[1]*nrm[1] + nrm[2]*nrm[2]);
		if(dv < FLT_EPSILON) dv = 1.;

		glBegin(GL_LINES);
			glVertex3d(vert[0], vert[1], vert[2]);
			glVertex3d(vert[0] + (size/dv)*nrm[0], vert[1] + (size/dv)*nrm[1], vert[2] + (size/dv)*nrm[2]);
		glEnd();
	}
}
/**-------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
void CHE_L0::draw_vf()
/*--------------------------------------------------------------------------*/
{
	/*
	float size = 0.1;
	float vert[3];
	float vf[3];

	glColor3d(0.5,0.5,0.5);		
	
	for (int i = 0 ; i < _nvert; i++){
		vert[0] = _G[i].x();
		vert[1] = _G[i].y();
		vert[2] = _G[i].z();

		 vf[0] = _G[i].fx();
		 vf[1] = _G[i].fy();
		 vf[2] = _G[i].fz();

		float dv = sqrt(vf[0]*vf[0] + vf[1]*vf[1] + vf[2]*vf[2]);
		if(dv < FLT_EPSILON) dv = 1.;

		glBegin(GL_LINES);
			glVertex3d(vert[0], vert[1], vert[2]);
			glVertex3d(vert[0] + (size/dv)*vf[0], vert[1] + (size/dv)*vf[1], vert[2] + (size/dv)*vf[2]);
		glEnd();
	}
	*/
}
/**-------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
void CHE_L0::read_ply( const char* file )
/*--------------------------------------------------------------------------*/
/** Reads a 3D triangulated model in the PLY file format.*/
{
  printf("CHE_L0::read_ply(%s)...", file) ;

  _G.clear();
  _V.clear();

  /*** the Ply object ***/
  PlyFile  *in_ply;
  PlyFace     face;
  PlyPoint       v;

  /** loop indices */
  int   j, ntrigs, nverts;
  int   i, elem_count;
  char  *elem_name;

  FILE *fp = fopen( file, "r" );

  if( fp==NULL ) { printf(" file not found.\n" ); return; }

  in_ply   = ::read_ply (fp);

  /** gets the number of faces and vertices*/
  ntrigs = nverts = 0 ;
  for ( i = 0; i < in_ply->num_elem_types; i++ )
  {
    elem_name = setup_element_read_ply ( in_ply, i, &elem_count );
    if ( equal_strings ( "vertex", elem_name ) )
      nverts = elem_count;
    if ( equal_strings ( "face",   elem_name ) )
      ntrigs = elem_count;
  }

  set_nvert(nverts);
  set_ntrig(ntrigs);

  _G.resize(  nverts  );
  _V.resize( 3*ntrigs );

  /* examine each element type that is in the file (PlyVertex, PlyFace) */
  for ( i = 0; i < in_ply->num_elem_types; i++ )
  {
    /* prepare to read the i'th list of elements */
    elem_name = setup_element_read_ply ( in_ply, i, &elem_count );

    if ( equal_strings ( "vertex", elem_name ) )
    {
      /* set up for getting PlyVertex elements */
      setup_property_ply ( in_ply, &plyvert_props[0] );
      setup_property_ply ( in_ply, &plyvert_props[1] );
      setup_property_ply ( in_ply, &plyvert_props[2] );
      setup_property_ply ( in_ply, &plyvert_props[3] );
      setup_property_ply ( in_ply, &plyvert_props[4] );
      setup_property_ply ( in_ply, &plyvert_props[5] );

      for ( j = 0; j < nverts; j++ )
      { 
		Vertex p;
        get_element_ply ( in_ply, ( void * ) &(v) );
		p.set_x ( (float)v.x  );
		p.set_y ( (float)v.y  );
		p.set_z ( (float)v.z  );
		p.set_nx( (float)v.nx );
		p.set_ny( (float)v.ny );
		p.set_nz( (float)v.nz );
		set_G(j, p);
      }
    }
    else if ( equal_strings ( "face", elem_name ) )
    {
      setup_property_ply ( in_ply, &plyface_props[0] );
      for ( j = 0; j < ntrigs; j++ )
      {
        get_element_ply ( in_ply, ( void * ) &face );
        if( face.nverts != 3 ) { printf("PLY importation: not a triangulated surface\n"); }
        _V[3*j]   = face.verts[0] ;
        _V[3*j+1] = face.verts[1] ;
        _V[3*j+2] = face.verts[2] ;
      }
    }
    else  /* all non-PlyVertex and non-PlyFace elements are grabbed here */
      get_other_element_ply ( in_ply );
  }
  close_ply( in_ply );
  free_ply ( in_ply );

  printf(" %d vertices and %d triangles found\n", nverts, ntrigs ) ;

  //Legalize the Model
  float min[3], max[3];

  compute_normals(); cout << endl;

  bounding_box  ( min, max);
 // legalize_model( min, max);
}
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
void CHE_L0::write_ply( const char* file, bool bin )
/*--------------------------------------------------------------------------*/
/** Writes a 3D triangulated model in the PLY file format.*/
{
  printf("CHE_L0::write_ply(%s)...", file) ;

  PlyFile    *ply;
  FILE       *fp = fopen( file, "w" );

  PlyFace     face ;
  int         v[3] ;
  char       *elem_names[]  = { "vertex", "face" };
  ply = ::write_ply ( fp, 2, elem_names, bin? PLY_BINARY_LE : PLY_ASCII );

  /* describe what properties go into the PlyVertex elements */
  describe_element_ply  ( ply, "vertex", nvert() );
  describe_property_ply ( ply, &plyvert_props[0] );
  describe_property_ply ( ply, &plyvert_props[1] );
  describe_property_ply ( ply, &plyvert_props[2] );
  describe_property_ply ( ply, &plyvert_props[3] );
  describe_property_ply ( ply, &plyvert_props[4] );
  describe_property_ply ( ply, &plyvert_props[5] );

  /* describe PlyFace properties (just list of PlyVertex indices) */
  describe_element_ply  ( ply, "face",   ntrig() );
  describe_property_ply ( ply, &plyface_props[0] );

  header_complete_ply ( ply );

  /* set up and write the PlyVertex elements */
  put_element_setup_ply ( ply, "vertex" );

  for( int i=0; i<nvert(); i++ )
  {
    Vertex &vt = _G[i] ;
	  PlyPoint p;

	  p.x = (float)vt.x();
	  p.y = (float)vt.y();
	  p.z = (float)vt.z();
	  p.nx = (float)vt.nx();
	  p.ny = (float)vt.ny();
	  p.nz = (float)vt.nz();

    put_element_ply ( ply, ( void * ) &(p) );
  }

  /* set up and write the PlyFace elements */
  put_element_setup_ply ( ply, "face" );
  face.nverts = 3 ;

  for( int i = 0 ; i < ntrig() ; ++i )
  {
    v[0] = _V[ 3*i ] ;
    v[1] = _V[3*i+1] ;
    v[2] = _V[3*i+2] ;
    face.verts = v ;
    put_element_ply ( ply, ( void * ) &face );
  }

  close_ply ( ply );
  free_ply ( ply );
  printf(" %d vertices and %d triangles written\n", nvert(), ntrig() ) ;
}
/*--------------------------------------------------------------------------*/



#include <cmath>
#include <Util/exceptions.h>
#include "scene.h"

using namespace Ray;
using namespace Util;
using namespace std;

///////////
// Scene //
///////////
Point3D Scene::Reflect( Point3D v , Point3D n )
{
	//////////////////
	// Reflect here //
	//////////////////
	WARN_ONCE( "method undefined" );
	return Point3D();
}

bool Scene::Refract( Point3D v , Point3D n , double ir , Point3D& refract )
{
	//////////////////
	// Refract here //
	//////////////////
	WARN_ONCE( "method undefined" );
	return false;
}

// called directly by the scene, which calls processFirstIntersection's
Point3D Scene::getColor( Ray3D ray , int rDepth , Point3D cLimit , unsigned int lightSamples , unsigned int tIdx )
{
	Point3D color;
	RayTracingStats::IncrementRayNum();
	ShapeProcessingInfo spInfo;
	RayIntersectionFilter rFilter = []( double ){ return true; };
	// sets color of intersected ray
	RayIntersectionKernel rKernel = [&]( const ShapeProcessingInfo &spInfo , const RayShapeIntersectionInfo &_iInfo )
	{
		/////////////////////////////////////////////////////////
		// Create the computational kernel that gets the color //
		/////////////////////////////////////////////////////////
		// WARN_ONCE( "method undefined" );
		color = Point3D( 1. , 1. , 1. ); // color white on hit
		return true;
	};

	processFirstIntersection( ray , BoundingBox1D( Epsilon , Infinity ) , rFilter , rKernel , spInfo , tIdx );

	// if ((ray.direction[0] - 0) < 0.1 && (ray.direction[1] + 0.707107) < 0.1 && (ray.direction[2] + 0.707107) < 0.1) {
	// 	cout << "ray in scene.todo: " << ray.direction << endl;
	// }

	return color;
}

//////////////
// Material //
//////////////
void Material::drawOpenGL( GLSLProgram * glslProgram ) const
{
	//////////////////////////////
	// Do OpenGL rendering here //
	//////////////////////////////
	WARN_ONCE( "method undefined" );

	// Sanity check to make sure that OpenGL state is good
	ASSERT_OPEN_GL_STATE();	
}

/////////////
// Texture //
/////////////
void Texture::initOpenGL( void )
{
	///////////////////////////////////
	// Do OpenGL texture set-up here //
	///////////////////////////////////
	WARN_ONCE( "method undefined" );

	// Sanity check to make sure that OpenGL state is good
	ASSERT_OPEN_GL_STATE();	
}

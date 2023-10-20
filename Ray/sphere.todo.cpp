#include <cmath>
#include <Util/exceptions.h>
#include "scene.h"
#include "sphere.h"

using namespace Ray;
using namespace Util;
using namespace std;

////////////
// Sphere //
////////////

void Sphere::init( const LocalSceneData &data )
{
	// Set the material pointer
	if( _materialIndex<0 ) THROW( "negative material index: " , _materialIndex );
	else if( _materialIndex>=data.materials.size() ) THROW( "material index out of bounds: " , _materialIndex , " <= " , data.materials.size() );
	else _material = &data.materials[ _materialIndex ];
	_primitiveNum = 1;

	///////////////////////////////////
	// Do any additional set-up here //
	///////////////////////////////////

	// get equation coefficients of sphere
	// P(x, y, z) = (x - a)^2 + (y - b)^2 + (z - c)^2 - r^2
	// P(x, y, z) = x^2 - 2ax + a^2 + y^2 - 2by + b^2 + z^2 - 2cz + c^2 - r^2
	// P(x, y, z) = x^2 - 2ax + y^2 - 2by + z^2 - 2cz + a^2 + b^2 + c^2 - r^2
	int a = center[0], b = center[1], c = center[2], r = radius;

	poly.coefficient(2u, 0u, 0u) = 1;
	poly.coefficient(1u, 0u, 0u) = -2 * a;
	poly.coefficient(0u, 2u, 0u) = 1;
	poly.coefficient(0u, 1u, 0u) = -2 * b;
	poly.coefficient(0u, 0u, 2u) = 1;
	poly.coefficient(0u, 0u, 1u) = -2 * c;
	poly.coefficient(0u, 0u, 0u) = pow(a, 2) + pow(b, 2) + pow(c, 2) - pow(r, 2);

}
void Sphere::updateBoundingBox( void )
{
	///////////////////////////////
	// Set the _bBox object here //
	///////////////////////////////
	WARN_ONCE( "method undefined" );
}
void Sphere::initOpenGL( void )
{
	///////////////////////////
	// Do OpenGL set-up here //
	///////////////////////////
	WARN_ONCE( "method undefined" );

	// Sanity check to make sure that OpenGL state is good
	ASSERT_OPEN_GL_STATE();	
}

// BoundingBox1D range is a bounding range of time — can't be too low for reflections, don't want to intersect from your starting point.
// RayIntersectionFilter is function: double -> bool
// RayIntersectionKernel is function: const ShapeProcessingInfo & , const RayShapeIntersectionInfo & -> void

// use a polynomial 3D, which you will make the equation for the sphere. You can run poly[ray] to get Polynomial1D intersection, which is the time of intersection I think.
bool Sphere::processFirstIntersection( const Ray3D &ray , const BoundingBox1D &range , const RayIntersectionFilter &rFilter , const RayIntersectionKernel &rKernel , ShapeProcessingInfo spInfo , unsigned int tIdx ) const
{
	RayTracingStats::IncrementRayPrimitiveIntersectionNum();
	spInfo.material = _material;

	//////////////////////////////////////////////////////////////
	// Compute the intersection of the sphere with the ray here //
	//////////////////////////////////////////////////////////////

	// get intersection with this sphere & given ray
	Util::Polynomial1D< 2 > intersection = poly(ray);
	double roots[2];
	unsigned int rootNum = intersection.roots(roots);

	// find closest intersection
	double left_bound = range[0][0], right_bound = range[1][0]; // first idx takes point<1> to double
	double closest_time = Infinity;
	for (int i = 0; i < rootNum; i++) {
		if (left_bound <= roots[i] && roots[i] <= right_bound && roots[i] < closest_time) {
			closest_time = roots[i];
		}
	}

	if (closest_time == Infinity) return false;

	// make intersection info and invoke rKernel
	RayShapeIntersectionInfo rsii = RayShapeIntersectionInfo();
	rsii.t = closest_time;
	rsii.position = ray(closest_time);
	Point3D normal = rsii.position - center;
	rsii.normal = normal / normal.length();
	rsii.texture = Point2D(0, 0); // TODO: figure out how to get texture coords
	rKernel(spInfo, rsii);
	return true;
}

int Sphere::processAllIntersections( const Ray3D &ray , const BoundingBox1D &range , const RayIntersectionFilter &rFilter , const RayIntersectionKernel &rKernel , ShapeProcessingInfo spInfo , unsigned int tIdx ) const
{
	RayTracingStats::IncrementRayPrimitiveIntersectionNum();
	spInfo.material = _material;

	//////////////////////////////////////////////////////////////
	// Compute the intersection of the sphere with the ray here //
	//////////////////////////////////////////////////////////////
	WARN_ONCE( "method undefined" );
	return 0;
}

bool Sphere::isInside( Point3D p ) const
{
	//////////////////////////////////////////////////////
	// Determine if the point is inside the sphere here //
	//////////////////////////////////////////////////////
	WARN_ONCE( "method undefined" );
	return false;
}

void Sphere::drawOpenGL( GLSLProgram * glslProgram ) const
{
	//////////////////////////////
	// Do OpenGL rendering here //
	//////////////////////////////
	WARN_ONCE( "method undefined" );

	// Sanity check to make sure that OpenGL state is good
	ASSERT_OPEN_GL_STATE();	
}

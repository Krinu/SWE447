/*
This whole project is a gigantic train wreck. It's cobbled together code from five or six different tutorials, using math that I
don't fully understand, often in ways that seem to me like they shouldn't work.

Well, that's C++ for you.

Going through the list of where I learned how to do various parts:
- R/T intersection, R/AABB intersection, PPM output, geometry.h: scratchapixel
- KD tree optimization: 

*/

#define NOMINMAX
#define M_PI 3.14159265358
#define EPSILON 1e-8
#define WIDTH 640
#define HEIGHT 480

#include "geometry.h"
#include <GL\glew.h>
#include <GLFW\glfw3.h>
#include <fstream>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <string>
#include <sstream>
#include <vector>
#include <time.h>

using namespace std;

inline float deg2rad(const float &deg) { return deg * M_PI / 180; }

inline float clamp(const float &lo, const float &hi, const float &v) { return max(lo, min(hi, v)); }

static Vec3f *fb = new Vec3f[HEIGHT * WIDTH]; // ugly but hey it works
static Vec3f *pixel = fb;
GLbyte imageData[HEIGHT * WIDTH * 3]; // color data for opengl
static int pixel_x = 0;
static int pixel_y = 0;

class aabb { // axis-aligned bounding box, for optimization
public:
	aabb();
	aabb(const Vec3f &b0, const Vec3f &b1) { bounds[0] = b0, bounds[1] = b1; }
	void expandBoundingBox(aabb boundingbox); // expands the called bounding box to hold the new one
	int getLongestAxis();
	Vec3f bounds[2];
};

aabb::aabb() {
	bounds[0] = 0; // min
	bounds[1] = 0; // max
}

void aabb::expandBoundingBox(aabb boundingbox) { // this may cause bugs due to sign checks on the box
	// std::cout << "bounding box expansion" << endl;
	// std::cout << "existing box: " << bounds[0].x << ", " << bounds[0].y << ", " << bounds[0].z << endl;
	// std::cout << "comparison box: " << boundingbox.bounds[0].x << ", " << boundingbox.bounds[0].y << ", " << boundingbox.bounds[0].z << endl;
	bounds[0].x = min(bounds[0].x, boundingbox.bounds[0].x);
	bounds[0].y = min(bounds[0].y, boundingbox.bounds[0].y);
	bounds[0].z = min(bounds[0].z, boundingbox.bounds[0].z);
	bounds[1].x = max(bounds[1].x, boundingbox.bounds[1].x);
	bounds[1].y = max(bounds[1].y, boundingbox.bounds[1].y);
	bounds[1].z = max(bounds[1].z, boundingbox.bounds[1].z);
	// std::cout << "new box: " << bounds[0].x << ", " << bounds[0].y << ", " << bounds[0].z << endl;
}

int aabb::getLongestAxis() {
	float xlength = abs(bounds[0].x) + abs(bounds[1].x);
	float ylength = abs(bounds[0].y) + abs(bounds[1].y);
	float zlength = abs(bounds[0].z) + abs(bounds[1].z);
	float biggest = max((xlength, ylength), zlength);

	if (xlength < biggest && ylength < biggest) {
		return 2; // x and y are smaller than biggest, Z is largest
	}
	else if (xlength < biggest && zlength < biggest) {
		return 1; // if x and z are smaller than biggest, Y is largest
	}
	else {
		return 0; // only other option is for X to be largest
	}

	// and since this function is only used for bounding box cutting for the KD tree, if it's a bit "fuzzy"
	// on the border case, it doesn't matter since it's really just a speed optimization for the renderer
	// anyway.
}

class triangle {
public:
	triangle(); // empty constructor
	aabb getBoundingBox();
	Vec3f getMidPoint();
	Vec3f v0 = { 0, 0, 0 };
	Vec3f v1 = { 0, 0, 0 };
	Vec3f v2 = { 0, 0, 0 };
	float t; // for holding texture/color coords
	float u;
	float v;
};

triangle::triangle() {
}

aabb triangle::getBoundingBox() { // MAKE ME A BOUNDING BOX
	// return aabb((min((v0.x, v1.x), v2.x), min((v0.y, v1.y), v2.y), min((v0.z, v1.z), v2.z)), (max((v0.x, v1.x), v2.x), max((v0.y, v1.y), v2.y), max((v0.z, v1.z), v2.z)));
	float smallestx;
	float smallesty;
	float smallestz;
	float largestx;
	float largesty;
	float largestz;

	// smallest
	smallestx = v0.x;
	if (v1.x < smallestx) {
		smallestx = v1.x;
	}
	if (v2.x < smallestx) {
		smallestx = v2.x;
	}

	smallesty = v0.y;
	if (v1.y < smallesty) {
		smallesty = v1.y;
	}
	if (v2.y < smallesty) {
		smallesty = v2.y;
	}

	smallestz = v0.z;
	if (v1.z < smallestz) {
		smallestz = v1.z;
	}
	if (v2.z < smallestz) {
		smallestz = v2.z;
	}

	// largest
	largestx = v0.x;
	if (v1.x > largestx) {
		largestx = v1.x;
	}
	if (v2.x > largestx) {
		largestx = v2.x;
	}

	largesty = v0.y;
	if (v1.y > largesty) {
		largesty = v1.y;
	}
	if (v2.y > largesty) {
		largesty = v2.y;
	}

	largestz = v0.z;
	if (v1.z > largestz) {
		largestz = v1.z;
	}
	if (v2.z > largestz) {
		largestz = v2.z;
	}

	Vec3f upper(largestx, largesty, largestz);
	Vec3f lower(smallestx, smallesty, smallestz);
	aabb box(lower, upper);
	return box;
}

Vec3f triangle::getMidPoint() {
	float x = (v0.x + v1.x + v2.x) / 3.0;
	float y = (v0.y + v1.y + v2.y) / 3.0;
	float z = (v0.z + v1.z + v2.z) / 3.0;
	return Vec3f(x, y, z);
}

class ray {
public:
	ray(const Vec3f &orig, const Vec3f &dir) : orig(orig), dir(dir)
	{
		invdir = 1 / dir;
		sign[0] = (invdir.x < 0);
		sign[1] = (invdir.y < 0);
		sign[2] = (invdir.z < 0);
	}
	void updateInvDir();
	Vec3f orig; // origin position
	Vec3f dir; // direction
	Vec3f invdir; // inverse direction
	int sign[3]; // sign of the components (used for bounding box checks)
};

void ray::updateInvDir() {
	invdir = 1 / dir;
	sign[0] = (invdir.x < 0);
	sign[1] = (invdir.y < 0);
	sign[2] = (invdir.z < 0);
}

void onHit(triangle target, ray r) {
	// *pixel = target.u * Vec3f(1, 0, 0) + target.v * Vec3f(0, 1, 0) + (1 - target.u - target.v) * Vec3f(0, 0, 1); // for color
	// *pixel++;
	pixel[pixel_x*pixel_y] = Vec3f(1, 1, 1); // white
}

bool doesIntersectBox(ray &r, aabb &box, float &t)
{
	r.updateInvDir(); // make sure our ray has its inverse determinant up to date

	float tmin, tmax, tymin, tymax, tzmin, tzmax;

	tmin = (box.bounds[r.sign[0]].x - r.orig.x) * r.invdir.x;
	tmax = (box.bounds[1 - r.sign[0]].x - r.orig.x) * r.invdir.x;
	tymin = (box.bounds[r.sign[1]].y - r.orig.y) * r.invdir.y;
	tymax = (box.bounds[1 - r.sign[1]].y - r.orig.y) * r.invdir.y;

	if ((tmin > tymax) || (tymin > tmax))
		return false;

	if (tymin > tmin)
		tmin = tymin;
	if (tymax < tmax)
		tmax = tymax;

	tzmin = (box.bounds[r.sign[2]].z - r.orig.z) * r.invdir.z;
	tzmax = (box.bounds[1 - r.sign[2]].z - r.orig.z) * r.invdir.z;

	if ((tmin > tzmax) || (tzmin > tmax))
		return false;

	if (tzmin > tmin)
		tmin = tzmin;
	if (tzmax < tmax)
		tmax = tzmax;

	t = tmin;

	if (t < 0) {
		t = tmax;
		if (t < 0) return false;
	}
	return true;
}

bool doesIntersectTriangle( // ray/triangle intersection
	const Vec3f &origin,
	const Vec3f &direction,
	const Vec3f &v0,
	const Vec3f &v1,
	const Vec3f &v2,
	float &t,
	float &u,
	float &v)
{
	Vec3f v0v1 = v1 - v0;
	Vec3f v0v2 = v2 - v0;
	Vec3f pvec = direction.crossProduct(v0v2);
	float determinant = v0v1.dotProduct(pvec);

	if (determinant < EPSILON) {
		return false;
	}

	float inversedeterminant = 1 / determinant;

	Vec3f tvec = origin - v0;
	u = tvec.dotProduct(pvec) * inversedeterminant;
	if (u < 0 || u > 1) {
		return false;
	}

	Vec3f qvec = tvec.crossProduct(v0v1);
	v = direction.dotProduct(qvec) * inversedeterminant;
	if (v < 0 || u + v > 1) {
		return false;
	}

	t = v0v2.dotProduct(qvec) * inversedeterminant;

	return true;
}

class kdnode {
public:
	kdnode* buildkdtree(vector<triangle> triangles, int depth);
	bool checkIntersection(kdnode* node, ray &r, float &t);
	kdnode* leftnode;
	kdnode* rightnode;
	aabb boundingbox;
	vector<triangle> triangles;
};

kdnode* kdnode::buildkdtree(vector<triangle> triangles, int depth) {
	/*
	for (int i = 0; i < triangles.size(); i++) {
		std::cout << triangles[i].v0.x << ", " << triangles[i].v0.y << ", " << triangles[i].v0.z << endl;
		std::cout << triangles[i].v1.x << ", " << triangles[i].v1.y << ", " << triangles[i].v1.z << endl;
		std::cout << triangles[i].v2.x << ", " << triangles[i].v2.y << ", " << triangles[i].v2.z << endl;
	}*/

	// declare the vars we'll need for moving down the array
	// std::cout << "entered top of buildkdtree" << endl;
	// std::cout << "triangle vector size: " << triangles.size() << ", depth " << depth << endl;
	kdnode* node = new kdnode();
	node->triangles = triangles;
	node->leftnode = NULL;
	node->rightnode = NULL;
	node->boundingbox = aabb();


	if (triangles.size() == 0) { // if there are no triangles here, return a pointer to the passed node
		// std::cout << "triangles size was zero, returned existing node." << endl;
		return node;


	}

	if (triangles.size() == 1) { // if there is only one triangle in this node
		// std::cout << "triangles size was one." << endl;
		node->boundingbox = triangles[0].getBoundingBox(); // set the bounding box for the node equal to the only triangle in the scene
		node->leftnode = new kdnode(); // create a new left and right node for use in a moment
		node->rightnode = new kdnode();
		node->leftnode->triangles = vector<triangle>(); // empty vectors to ensure tree consistency
		node->rightnode->triangles = vector<triangle>();
		return node;
	}

	// if we've reached this point, we have at least two items left to add to our tree, which means we have work to do
	// first, set up the bounding box for the triangles - start with one, and then grow it
	node->boundingbox = triangles[0].getBoundingBox();

	// now grow the box to cover all triangles in the vector - counter starts at 1, since we did the 0th in the above line
	for (int i = 1; i < triangles.size(); i++) {
		node->boundingbox.expandBoundingBox(triangles[i].getBoundingBox());
	}

	// the bounding box now covers all the triangles in the vector. the next step is to figure out where to "cut" the structure
	// in order to split the tree. there are several ways to do this, with different advantages and disadvantages. Here, I just
	// split at the midpoint of the triangle group for ease of implementation.
	Vec3f midpoint(0, 0, 0);

	// now that we have a zero point listed, we add up each triangle's contribution to the average to get the midpoint.
	for (int i = 0; i < triangles.size(); i++) {
		midpoint = midpoint + (triangles[i].getMidPoint() * (1.0 / triangles.size()));
	}

	// with the midpoint in hand, we can split the tree around the longest axis (which is marginally better than X->Y->Z)
	vector<triangle> left_triangles;
	vector<triangle> right_triangles;

	int axis = node->boundingbox.getLongestAxis();
	for (int i = 0; i < triangles.size(); i++) {
		switch (axis) {
		case 0:
			midpoint.x >= triangles[i].getMidPoint().x ? right_triangles.push_back(triangles[i]) : left_triangles.push_back(triangles[i]);
			break;
		case 1:
			midpoint.y >= triangles[i].getMidPoint().y ? right_triangles.push_back(triangles[i]) : left_triangles.push_back(triangles[i]);
			break;
		case 2:
			midpoint.z >= triangles[i].getMidPoint().z ? right_triangles.push_back(triangles[i]) : left_triangles.push_back(triangles[i]);
			break;
		}
	}

	// we now have our split point. from here, we have to decide how far to recurse - some recommend stopping after distributing
	// a certain percentage of the triangles, others after a finite number of splits. For now, I'm going with 20 splits. After
	// hitting 20, it dumps the remaining triangles in the last left and right nodes and calls it a day. Note that it is possible
	// for this function to end earlier (see the triangle==0 and triangle==1 section up above), thus never triggering the max
	// depth end condition.
	
	if (depth < 3) { // TODO: DEBUG ME!!!
		node->leftnode = buildkdtree(left_triangles, depth + 1);
		node->rightnode = buildkdtree(right_triangles, depth + 1);
	}



	else {
		node->leftnode = new kdnode();
		node->rightnode = new kdnode();
		node->leftnode->triangles = vector<triangle>();
		node->rightnode->triangles = vector<triangle>();
	}

	return node;

}

bool kdnode::checkIntersection(kdnode* node, ray &r, float &t) {
	/*
	std::cout << "checkIntersection debug code" << endl;
	std::cout << "node triangles array size: " << node->triangles.size() << endl;
	std::cout << "ray origin: " << r.orig.x << ", " << r.orig.y << ", " << r.orig.z << endl;
	std::cout << "ray direction: " << r.dir.x << ", " << r.dir.y << ", " << r.dir.z << endl;
	std::cout << "contents of t: " << t << endl;
	std::cout << "bound box lower target: " << node->boundingbox.bounds[0].x << ", " << node->boundingbox.bounds[0].y << ", " << node->boundingbox.bounds[0].z << endl;
	std::cout << "bound box upper target: " << node->boundingbox.bounds[1].x << ", " << node->boundingbox.bounds[1].y << ", " << node->boundingbox.bounds[1].z << endl;
	*/
	if (doesIntersectBox(r, node->boundingbox, t)) {
		//std::cout << "left size: " << node->leftnode->triangles.size() << endl;
		//std::cout << "right size: " << node->rightnode->triangles.size() << endl;
		//std::cout << "ray hit bounding box" << endl;
		bool hit_triangle = false;


		if (node->leftnode->triangles.size() > 0 || node->rightnode->triangles.size() > 0) {
			// std::cout << "left or right child was greater than zero" << endl;
			bool hitleft = checkIntersection(node->leftnode, r, t);
			bool hitright = checkIntersection(node->rightnode, r, t);
			return hitleft || hitright;
		}
		else {
			// std::cout << "left or right child was NOT greater than zero" << endl;
			Vec3f collision();
			for (int i = 0; i < node->triangles.size(); i++) {
				if (doesIntersectTriangle(r.orig, r.dir, node->triangles[i].v0, node->triangles[i].v1, node->triangles[i].v2, node->triangles[i].t, node->triangles[i].u, node->triangles[i].v)) {
					//onHit(triangles[i], r);
					hit_triangle = true;
				}
			}
			if (hit_triangle) {
				return true;
			}
			return false;
		}
	}
	return false;
}



void loadScene(const char* filename, vector<triangle> &triangles) { 
	// adapted from https://en.wikibooks.org/wiki/OpenGL_Programming/Modern_OpenGL_Tutorial_Load_OBJ
	
	// DEBUG
	std::cout << "Loading model file " << filename << " from disk.\n";
	
	vector<Vec3f> vertices;
	vector<int> elements;
	int vertextracker = 0; // technically, I could use the modulo function here but meh

	// standard file i/o stuff
	ifstream in(filename, ios::in);
	
	if (!in)
	{
		cerr << "Cannot open " << filename << endl; exit(1);
	}


	string line;
	while (getline(in, line)) {
		if (line.substr(0, 2) == "v ") { 
			// build the vertex index. see comment below for what we use this for.
			istringstream s(line.substr(2));
			Vec3f v;
			s >> v.x;
			s >> v.y;
			s >> v.z;
			vertices.push_back(v);
			cout << v.x << ", " << v.y << ", " << v.z << endl;
		}
		else if (line.substr(0, 2) == "f ") {
			// faces in obj files are in the form of a line header ("f " in this case) followed by three values separated by spaces
			// while the face entry can have a lot of different pieces of data (with or without texure coordinate indices, etc.)
			// it doesn't change the fact that the first entry after "f " and before a second space is the first vertex
			// so while we might have "f 123//52" or "f 123/22" and those second values mean different things, it doesn't
			// change the fact that the vertex for that face is at index 123. We thus want to take the first ints we find and stop
			// reading at the first forward slash.

			// TODO: Replace with something that DOESN'T use three different kinds of string data.
			istringstream s(line.substr(2));
			string foo = s.str();
			cout << foo << endl;
			const char* tmpcstring = foo.c_str();
			int a;
			int b;
			int c;
			int unused;
			sscanf_s(tmpcstring, "%d//%d %d//%d %d//%d", &a, &unused, &b, &unused, &c, &unused);
			s >> a;
			s >> b;
			s >> c;
			a--;
			b--;
			c--;
			elements.push_back(a);
			elements.push_back(b);
			elements.push_back(c);
			cout << a << ", " << b << ", " << c << endl;
		}
	}


	// DEBUG
	std::cout << "Finished loading model file to vertex array.\n";

	// now that we have vertex data, turn it into triangles
// there are likely more elegant ways to do this, but this works
triangle t;
for (int i = 0; i < elements.size(); i++) {
	if (vertextracker == 0) {
		t.v0 = vertices[elements[i]];
		vertextracker++;
	}
	else if (vertextracker == 1) {
		t.v1 = vertices[elements[i]];
		vertextracker++;
	}
	else if (vertextracker == 2) {
		t.v2 = vertices[elements[i]];
		vertextracker = 0;
		triangles.push_back(t);
	}
}

// DEBUG
std::cout << "Finished converting vertex vector to triangle vector.\n";

}

void moveScene(float xdir, float ydir, float zdir, vector<triangle> &scene) {
	// DEBUG
	std::cout << "Relocating scene.\n";
	for (int i = 0; i < scene.size(); i++) {
		scene[i].v0.x += xdir;
		scene[i].v0.y += ydir;
		scene[i].v0.z += zdir;
		scene[i].v1.x += xdir;
		scene[i].v1.y += ydir;
		scene[i].v1.z += zdir;
		scene[i].v2.x += xdir;
		scene[i].v2.y += ydir;
		scene[i].v2.z += zdir;
	}
}

void scaleScene(float xscale, float yscale, float zscale, vector<triangle> &scene) {
	for (int i = 0; i < scene.size(); i++) {
		scene[i].v0.x *= xscale;
		scene[i].v0.y *= yscale;
		scene[i].v0.z *= zscale;
		scene[i].v1.x *= xscale;
		scene[i].v1.y *= yscale;
		scene[i].v1.z *= zscale;
		scene[i].v2.x *= xscale;
		scene[i].v2.y *= yscale;
		scene[i].v2.z *= zscale;
	}
}

void printScene(vector<triangle> &scene) {
	for (int i = 0; i < scene.size(); i++) {
		std::cout << "TRIANGLE " << (i + 1) << ":" << endl;
		std::cout << "   v0: " << scene[i].v0.x << ", " << scene[i].v0.y << ", " << scene[i].v0.z << endl;
		std::cout << "   v1: " << scene[i].v1.x << ", " << scene[i].v1.y << ", " << scene[i].v1.z << endl;
		std::cout << "   v2: " << scene[i].v2.x << ", " << scene[i].v2.y << ", " << scene[i].v2.z << endl;
	}
}

void tracescene_old() {
	Vec3f cols[3] = { { float(1.0), float(0.0), float(0.0) },{ float(0.0), float(1.0), float(0.0) },{ float(0.0), float(0.0), float(1.0) } };
	Vec3f *framebuffer = new Vec3f[WIDTH * HEIGHT];
	Vec3f *pix = framebuffer;

	float fov = 51.52;
	float scale = tan(deg2rad(fov * 0.5));
	float imageAspectRatio = WIDTH / (float)HEIGHT;

	vector<triangle> tris;
	loadScene("test.obj", tris);
	moveScene(0, 0, -5, tris);
	//scaleScene(1, -1, 1, tris);
	printScene(tris);

	Vec3f orig(0);

	int totalpixels = HEIGHT * WIDTH;
	int framesdrawn = 0;

	for (uint32_t j = 0; j < HEIGHT; ++j) {
		std::cout << "\b\b\b\b\b\b" << framesdrawn;
		for (uint32_t i = 0; i < WIDTH; ++i) {
			// DEBUG
			// compute primary ray
			float x = (2 * (i + 0.5) / (float)WIDTH - 1) * imageAspectRatio * scale;
			float y = (1 - 2 * (j + 0.5) / (float)HEIGHT) * scale;
			Vec3f dir(x, y, -1);
			dir.normalize();
			ray r(Vec3f(0, 0, 0), dir);

			for (int i = 0; i < tris.size(); i++) {
				if (doesIntersectTriangle(orig, dir, tris[i].v0, tris[i].v1, tris[i].v2, tris[i].t, tris[i].u, tris[i].v)) {
					// the following is SUPER UGLY DEPTH-BASED FLAT SHADING CODE
					// which takes advantage of the fact that we store color in
					// the same Vec3f format we use for vertices before writing
					// it out to the framebuffer.

					// Take the midpoint of the triangle we just hit
					Vec3f center = tris[i].getMidPoint();

					// Set the X and Y values to 1 and 0 (for coloring). If we
					// want to get fancy, we can swap vertices around to let us
					// make different colors.
					// center.x = 1;
					// center.y = 0;
					center.z = 1;

					// Normalize it so we have a spread from 0->1 for our depth
					// value (to prevent oversaturation, which would just get
					// clamped by loop that writes the color values).
					*pix = center.normalize();
					


				}
			}
			pix++;
			framesdrawn++;
		}
	}

	std::ofstream ofs("./out.ppm", std::ios::out | std::ios::binary);
	ofs << "P6\n" << WIDTH << " " << HEIGHT << "\n255\n";
	for (uint32_t i = 0; i < HEIGHT * WIDTH; ++i) {
		char r = (char)(255 * clamp(0, 1, framebuffer[i].x));
		char g = (char)(255 * clamp(0, 1, framebuffer[i].y));
		char b = (char)(255 * clamp(0, 1, framebuffer[i].z));
		imageData[(i * 3)] = r;
		imageData[(i * 3) + 1] = g;
		imageData[(i * 3) + 2] = b;
		//imageData[(i * 4) + 3] = 255;
		
		ofs << r << g << b;
	}

	ofs.close();
}

void printkdtree(kdnode* tree, int depth) {
	std::cout << "tree triangle size at depth " << depth << ": ";
	std::cout << tree->triangles.size() << endl;
	if (depth < 20) {
		if (tree->leftnode) {
			printkdtree(tree->leftnode, depth + 1);
		}
		if (tree->rightnode) {
			printkdtree(tree->rightnode, depth + 1);
		}
	}
}

void tracescene_new() {
	Vec3f cols[3] = { { float(1.0), float(0.0), float(0.0) },{ float(0.0), float(1.0), float(0.0) },{ float(0.0), float(0.0), float(1.0) } };
	Vec3f *framebuffer = new Vec3f[WIDTH * HEIGHT];
	Vec3f *pix = framebuffer;

	//std::cout << "debug..." << endl;

	for (int i = 0; i < (HEIGHT * WIDTH); i++) {
		pixel = new Vec3f(1, 0, 0);
	}

	//std::cout << "debug!!!!!!!!!!!" << endl;

	aabb testbox(Vec3f(-1, -1, -4), Vec3f(1, 1, -5));

	float fov = 51.52;
	float scale = tan(deg2rad(fov * 0.5));
	float imageAspectRatio = WIDTH / (float)HEIGHT;

	Vec3f origin(0);
	int hitcount = 0;

	// TEMP TESTING CODE
	vector<triangle> scene; // vector of tris for scene
	triangle testtri;
	loadScene("test.obj", scene); // initialize the scene
	moveScene(0, 0, -5, scene); // move the scene away from the origin
	kdnode* nodeptr = new kdnode();
	nodeptr = nodeptr->buildkdtree(scene, 0);
	std::cout << "scene built!" << endl;
	printkdtree(nodeptr, 0);
	int framesdrawn = 0;
	aabb testbox2(Vec3f(-1, -1, -4), Vec3f(1, 1, -5));

	//printkdtree(nodeptr, 0);

	for (int i = 0; i < HEIGHT; ++i) {
		std::cout << "\b\b\b\b\b\b" << framesdrawn;
		pixel_y = i;
		for (int j = 0; j < WIDTH; ++j) {
			pixel_x = j;
			float x = (2 * (j + 0.5) / (float)WIDTH - 1) * imageAspectRatio * scale; // x component of the primary ray, scaled for screen size
			float y = (1 - 2 * (i + 0.5) / (float)HEIGHT) * scale; // y component of the primary ray, also scaled for screen size
			Vec3f direction(x, y, -1); // where we aim the primary ray
			direction.normalize(); // normalize the vector

			ray primary(Vec3f(0, 0, 0), direction); // create primary ray
			float t; // I forget why I used this, I think it was for logging.
			/*
			if (doesIntersectBox(primary, testbox2, t)) {
				onHit(testtri, primary);
				hitcount++;
				// std::cout << "hit at " << pixel_x << ", " << pixel_y << endl;
			}*/
			
			if (nodeptr->checkIntersection(nodeptr, primary, t)) {
				hitcount++;
				// std::cout << "hit" << endl;
			}

			framesdrawn++;
		}
	}

	/*
	for (int i = 0; i < intersectedtris.size(); i++) {
		*pix = intersectedtris[i]->u * cols[0] + intersectedtris[i]->v * cols[1] + (1 - intersectedtris[i]->u - intersectedtris[i]->v) * cols[2];
		pix++;

	}*/
	/*
	std::ofstream ofs("./out_new.ppm", std::ios::out | std::ios::binary);
	ofs << "P6\n" << WIDTH << " " << HEIGHT << "\n255\n";
	for (uint32_t i = 0; i < HEIGHT * WIDTH; ++i) {
		char r = (char)(255 * clamp(0, 1, pixel[i].x));
		char g = (char)(255 * clamp(0, 1, pixel[i].y));
		char b = (char)(255 * clamp(0, 1, pixel[i].z));
		ofs << r << g << b;
	}

	ofs.close();
	*/
	/*
	for (int i = 0; i < HEIGHT; ++i) {
		for (int j = 0; j < WIDTH; ++j) {
			float x = (2 * (j + 0.5) / (float)WIDTH - 1) * imageAspectRatio * scale;
			float y = (1 - 2 * (i + 0.5) / (float)HEIGHT) * scale;
			Vec3f direction(x, y, -1);
			direction.normalize();

			float t;
			ray r(Vec3f(0), direction);


			if (doesIntersectBox(r, testbox, t)) {
				Vec3f hitloc = r.orig + r.dir * t;
				std::cout << r.orig << " " << hitloc << endl;
				hitcount++;
			}
		}
	}
	*/

	std::cout << endl << "out of " << (HEIGHT * WIDTH) << " rays, " << hitcount << " hit the box.\n";
	
}



int main()
{
	
	tracescene_new();
	tracescene_old();

	if (!glfwInit()) {
		fprintf(stderr, "ERROR: Failed to initialize GLFW.\n");
	}

	GLFWwindow* window;

	window = glfwCreateWindow(WIDTH, HEIGHT, "Raytracer 3.0: Visual Studio Edition", NULL, NULL);

	if (window == NULL) {
		fprintf(stderr, "ERROR: Failed to create GLFW window.\n");
		glfwTerminate();
		return -1;
	}

	glfwMakeContextCurrent(window); // make it so, captain
	glewExperimental = true;
	
	if (glewInit() != GLEW_OK) {
		fprintf(stderr, "ERROR: Failed to initialize GLEW.\n");
		return -1;
	}

	glfwSetInputMode(window, GLFW_STICKY_KEYS, GL_TRUE);

	// BEGIN ANNOYING SETUP
	// build and bind the texture from our image data
	GLuint texture;
	glGenTextures(1, &texture);
	glBindTexture(GL_TEXTURE_2D, texture);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, WIDTH, HEIGHT, 0, GL_RGB, GL_UNSIGNED_BYTE, imageData);
	glGenerateMipmap(GL_TEXTURE_2D);
	glBindTexture(GL_TEXTURE_2D, 0);

	// assign the texture to a FBO (which is actually a reference to the existing framebuffer)
	// thank you stackoverflow (http://stackoverflow.com/questions/30488155/opengl-fastest-way-to-draw-2d-image)
	GLuint readFboId;
	glGenFramebuffers(1, &readFboId);
	glBindFramebuffer(GL_READ_FRAMEBUFFER, readFboId);
	glFramebufferTexture2D(GL_READ_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, texture, 0);
	glBindFramebuffer(GL_READ_FRAMEBUFFER, 0);

	
	do {
		glfwPollEvents();

		// and here the magic happens - well, as magic as this can be
		// draw it out
		glBindFramebuffer(GL_READ_FRAMEBUFFER, readFboId);
		glBlitFramebuffer(0, 0, WIDTH, HEIGHT, 0, 0, WIDTH, HEIGHT, GL_COLOR_BUFFER_BIT, GL_LINEAR);
		glBindFramebuffer(GL_READ_FRAMEBUFFER, 0);

		glfwSwapBuffers(window);
	} while (glfwGetKey(window, GLFW_KEY_ESCAPE) != GLFW_PRESS && glfwWindowShouldClose(window) == 0);


}
#include <stdio.h>
#include <stdlib.h>

#include "ml6.h"
#include "display.h"
#include "draw.h"
#include "matrix.h"
#include "math.h"
#include "gmath.h"


void swap(double * a, double * b){
	double t = *a;
	*a = *b;
	*b = t;
}


/*======== void scanline_convert() ==========
   Inputs: struct matrix *points
          int i
          screen s
          zbuffer zb
   Returns:

   Fills in polygon i by drawing consecutive horizontal (or vertical) lines.

   Color should be set differently for each polygon.
   ====================*/
void scanline_convert( struct matrix *points, int i, screen s, zbuffer zbuf ) {

	color c;
	c.red = rand() % 256;
	c.green = rand() % 256;
	c.blue = rand() % 256;

	double xb = points->m[0][i];
	double yb = points->m[1][i];
	double zb = points->m[2][i];

	double xm = points->m[0][i+1];
	double ym = points->m[1][i+1];
	double zm = points->m[2][i+1];

	double xt = points->m[0][i+2];
	double yt = points->m[1][i+2];
	double zt = points->m[2][i+2];

    printf("%f %f %f\n", yb, ym, yt);

    if (yb > ym) {
		swap(&xb, &xm);
		swap(&yb, &ym);
		swap(&zb, &zm);
	}

	if (ym > yt) {
		swap(&xm, &xt);
		swap(&ym, &yt);
		swap(&zm, &zt);
	}

    if (yb > ym) {
		swap(&xb, &xm);
		swap(&yb, &ym);
		swap(&zb, &zm);
	}

    printf("%lf %lf %lf\n", yb, ym, yt);

	double x0 = xb, x1 = xb;
	double z0 = zb, z1 = zb;
	double y = yb;

	double dx0 = (xt - xb) / (yt - yb + 1);
	double dx1_0 = (xm - xb) / (ym - yb + 1);
	double dx1_1 = (xt - xm) / (yt - ym + 1);

	double dz0 = (zt - zb) / (yt - yb + 1);
	double dz1_0 = (zm - zb) / (ym - yb + 1);
	double dz1_1 = (zt - zm) / (yt - ym + 1);

	double dx1 = dx1_0;
	double dz1 = dz1_0;

	if((int) ym == (int) yb) {
		dx1 = dx1_1;
		x1 = xm;

		dz1 = dz1_1;
		z1 = zm;
	}

    // printf("x: %lf %lf\n", x0, x1);
    // printf("y: %lf\n", y);
    // printf("z: %lf %lf\n", z0, z1);
    //
    // printf("dx: %lf %lf %lf\n", dx0, dx1_0, dx1_1);
    // printf("dz: %lf %lf %lf\n", dz0, dz1_0, dz1_1);

	while(y <= yt) {
		draw_line(x0, y, z0, x1, y, z1, s, zbuf, c);

		x0 += dx0;
		x1 += dx1;

		z0 += dz0;
		z1 += dz1;

		y++;

		if ((int) y == (int) ym) {
            printf("switch\n");
			dx1 = dx1_1;
			x1 = xm;

			dz1 = dz1_1;
			z1 = zm;
		}

	}




}

/*======== void add_polygon() ==========
   Inputs:   struct matrix *polygons
            double x0
            double y0
            double z0
            double x1
            double y1
            double z1
            double x2
            double y2
            double z2
   Returns:
   Adds the vertices (x0, y0, z0), (x1, y1, z1)
   and (x2, y2, z2) to the polygon matrix. They
   define a single triangle surface.
   ====================*/
void add_polygon( struct matrix *polygons,
                  double x0, double y0, double z0,
                  double x1, double y1, double z1,
                  double x2, double y2, double z2 ) {
	add_point(polygons, x0, y0, z0);
	add_point(polygons, x1, y1, z1);
	add_point(polygons, x2, y2, z2);
}

/*======== void draw_polygons() ==========
   Inputs:   struct matrix *polygons
            screen s
            color c
   Returns:
   Goes through polygons 3 points at a time, drawing
   lines connecting each points to create bounding triangles
   ====================*/
void draw_polygons( struct matrix *polygons, screen s, zbuffer zb, color c ) {
	if ( polygons->lastcol < 3 ) {
		printf("Need at least 3 points to draw a polygon!\n");
		return;
	}

	int point;
	double *normal;

	for (point=0; point < polygons->lastcol-2; point+=3) {

		normal = calculate_normal(polygons, point);

		if ( normal[2] > 0 ) {

            printf("filling in triangle %d\n", point);
			scanline_convert(polygons, point, s, zb);

			// draw_line( polygons->m[0][point],
			//            polygons->m[1][point],
			//            polygons->m[2][point],
			//            polygons->m[0][point+1],
			//            polygons->m[1][point+1],
			//            polygons->m[2][point+1],
			//            s, zb, c);
			// draw_line( polygons->m[0][point+2],
			//            polygons->m[1][point+2],
			//            polygons->m[2][point+2],
			//            polygons->m[0][point+1],
			//            polygons->m[1][point+1],
			//            polygons->m[2][point+1],
			//            s, zb, c);
            //
			// draw_line( polygons->m[0][point],
			//            polygons->m[1][point],
			//            polygons->m[2][point],
			//            polygons->m[0][point+2],
			//            polygons->m[1][point+2],
			//            polygons->m[2][point+2],
			//            s, zb, c);


		}
	}
}

/*======== void add_box() ==========
   Inputs:   struct matrix * edges
            double x
            double y
            double z
            double width
            double height
            double depth

   add the points for a rectagular prism whose
   upper-left-front corner is (x, y, z) with width,
   height and depth dimensions.
   ====================*/
void add_box( struct matrix *polygons,
              double x, double y, double z,
              double width, double height, double depth ) {
	double x0, y0, z0, x1, y1, z1;
	x0 = x;
	x1 = x+width;
	y0 = y;
	y1 = y-height;
	z0 = z;
	z1 = z-depth;


	//front
	add_polygon(polygons, x, y, z, x1, y1, z, x1, y, z);
	add_polygon(polygons, x, y, z, x, y1, z, x1, y1, z);
	//back
	add_polygon(polygons, x1, y, z1, x, y1, z1, x, y, z1);
	add_polygon(polygons, x1, y, z1, x1, y1, z1, x, y1, z1);

	//right side
	add_polygon(polygons, x1, y, z, x1, y1, z1, x1, y, z1);
	add_polygon(polygons, x1, y, z, x1, y1, z, x1, y1, z1);
	//left side
	add_polygon(polygons, x, y, z1, x, y1, z, x, y, z);
	add_polygon(polygons, x, y, z1, x, y1, z1, x, y1, z);

	//top
	add_polygon(polygons, x, y, z1, x1, y, z, x1, y, z1);
	add_polygon(polygons, x, y, z1, x, y, z, x1, y, z);
	//bottom
	add_polygon(polygons, x, y1, z, x1, y1, z1, x1, y1, z);
	add_polygon(polygons, x, y1, z, x, y1, z1, x1, y1, z1);
}


/*======== void add_sphere() ==========
   Inputs:   struct matrix * points
            double cx
            double cy
            double cz
            double r
            int step

   adds all the points for a sphere with center (cx, cy, cz)
   and radius r using step points per circle/semicircle.

   Since edges are drawn using 2 points, add each point twice,
   or add each point and then another point 1 pixel away.

   should call generate_sphere to create the necessary points
   ====================*/
void add_sphere( struct matrix * edges,
                 double cx, double cy, double cz,
                 double r, int step ) {

	struct matrix *points = generate_sphere(cx, cy, cz, r, step);
	int p0, p1, p2, p3, lat, longt;
	int latStop, longStop, latStart, longStart;
	latStart = 0;
	latStop = step;
	longStart = 1;
	longStop = step;

	//step++; needed for my triangles
	for ( lat = latStart; lat < latStop; lat++ ) {
		for ( longt = longStart; longt < longStop; longt++ ) {

			/*Milan's Triangles*/
			p0 = lat * (step+1) + longt;
			p1 = p0 + 1;
			p2 = (p1 + step) % (step * (step+1));
			p3 = (p0 + step) % (step * (step+1));

			add_polygon( edges, points->m[0][p0],
			             points->m[1][p0],
			             points->m[2][p0],
			             points->m[0][p1],
			             points->m[1][p1],
			             points->m[2][p1],
			             points->m[0][p2],
			             points->m[1][p2],
			             points->m[2][p2]);
			add_polygon( edges, points->m[0][p0],
			             points->m[1][p0],
			             points->m[2][p0],
			             points->m[0][p2],
			             points->m[1][p2],
			             points->m[2][p2],
			             points->m[0][p3],
			             points->m[1][p3],
			             points->m[2][p3]);

			/*My Triangles*/
			/* p0 = lat * (step) + longt; */
			/* p1 = p0+1; */
			/* p2 = (p1+step) % (step * (step-1)); */
			/* p3 = (p0+step) % (step * (step-1)); */

			/* //printf("p0: %d\tp1: %d\tp2: %d\tp3: %d\n", p0, p1, p2, p3); */
			/* if (longt < step - 2) */
			/*   add_polygon( edges, points->m[0][p0], */
			/*                points->m[1][p0], */
			/*                points->m[2][p0], */
			/*                points->m[0][p1], */
			/*                points->m[1][p1], */
			/*                points->m[2][p1], */
			/*                points->m[0][p2], */
			/*                points->m[1][p2], */
			/*                points->m[2][p2]); */
			/* if (longt > 0 ) */
			/*   add_polygon( edges, points->m[0][p0], */
			/*                points->m[1][p0], */
			/*                points->m[2][p0], */
			/*                points->m[0][p2], */
			/*                points->m[1][p2], */
			/*                points->m[2][p2], */
			/*                points->m[0][p3], */
			/*                points->m[1][p3], */
			/*                points->m[2][p3]); */
		}
	}
	free_matrix(points);
}

/*======== void generate_sphere() ==========
   Inputs:   struct matrix * points
            double cx
            double cy
            double cz
            double r
            int step
   Returns: Generates all the points along the surface
           of a sphere with center (cx, cy, cz) and
           radius r using step points per circle/semicircle.
           Returns a matrix of those points
   ====================*/
struct matrix * generate_sphere(double cx, double cy, double cz,
                                double r, int step ) {

	struct matrix *points = new_matrix(4, step * step);
	int circle, rotation, rot_start, rot_stop, circ_start, circ_stop;
	double x, y, z, rot, circ;

	rot_start = 0;
	rot_stop = step;
	circ_start = 0;
	circ_stop = step;

	for (rotation = rot_start; rotation < rot_stop; rotation++) {
		rot = (double)rotation / step;

		for(circle = circ_start; circle <= circ_stop; circle++) {
			circ = (double)circle / step;

			x = r * cos(M_PI * circ) + cx;
			y = r * sin(M_PI * circ) *
			    cos(2*M_PI * rot) + cy;
			z = r * sin(M_PI * circ) *
			    sin(2*M_PI * rot) + cz;

			/* printf("rotation: %d\tcircle: %d\n", rotation, circle); */
			/* printf("rot: %lf\tcirc: %lf\n", rot, circ); */
			/* printf("sphere point: (%0.2f, %0.2f, %0.2f)\n\n", x, y, z); */
			add_point(points, x, y, z);
		}
	}

	return points;
}

/*======== void add_torus() ==========
   Inputs:   struct matrix * points
            double cx
            double cy
            double cz
            double r1
            double r2
            double step
   Returns:

   adds all the points required for a torus with center (cx, cy, cz),
   circle radius r1 and torus radius r2 using step points per circle.

   should call generate_torus to create the necessary points
   ====================*/
void add_torus( struct matrix * edges,
                double cx, double cy, double cz,
                double r1, double r2, int step ) {

	struct matrix *points = generate_torus(cx, cy, cz, r1, r2, step);
	int p0, p1, p2, p3, lat, longt;
	int latStop, longStop, latStart, longStart;
	latStart = 0;
	latStop = step;
	longStart = 0;
	longStop = step;

	for ( lat = latStart; lat < latStop; lat++ ) {
		for ( longt = longStart; longt < longStop; longt++ ) {
			p0 = lat * step + longt;
			if (longt == step - 1)
				p1 = p0 - longt;
			else
				p1 = p0 + 1;
			p2 = (p1 + step) % (step * step);
			p3 = (p0 + step) % (step * step);

			//printf("p0: %d\tp1: %d\tp2: %d\tp3: %d\n", p0, p1, p2, p3);
			add_polygon( edges, points->m[0][p0],
			             points->m[1][p0],
			             points->m[2][p0],
			             points->m[0][p3],
			             points->m[1][p3],
			             points->m[2][p3],
			             points->m[0][p2],
			             points->m[1][p2],
			             points->m[2][p2]);
			add_polygon( edges, points->m[0][p0],
			             points->m[1][p0],
			             points->m[2][p0],
			             points->m[0][p2],
			             points->m[1][p2],
			             points->m[2][p2],
			             points->m[0][p1],
			             points->m[1][p1],
			             points->m[2][p1]);
		}
	}
	free_matrix(points);
}

/*======== void generate_torus() ==========
   Inputs:   struct matrix * points
            double cx
            double cy
            double cz
            double r
            int step
   Returns: Generates all the points along the surface
           of a torus with center (cx, cy, cz),
           circle radius r1 and torus radius r2 using
           step points per circle.
           Returns a matrix of those points
   ====================*/
struct matrix * generate_torus( double cx, double cy, double cz,
                                double r1, double r2, int step ) {

	struct matrix *points = new_matrix(4, step * step);
	int circle, rotation, rot_start, rot_stop, circ_start, circ_stop;
	double x, y, z, rot, circ;

	rot_start = 0;
	rot_stop = step;
	circ_start = 0;
	circ_stop = step;

	for (rotation = rot_start; rotation < rot_stop; rotation++) {
		rot = (double)rotation / step;

		for(circle = circ_start; circle < circ_stop; circle++) {
			circ = (double)circle / step;

			x = cos(2*M_PI * rot) *
			    (r1 * cos(2*M_PI * circ) + r2) + cx;
			y = r1 * sin(2*M_PI * circ) + cy;
			z = -1*sin(2*M_PI * rot) *
			    (r1 * cos(2*M_PI * circ) + r2) + cz;

			//printf("rotation: %d\tcircle: %d\n", rotation, circle);
			//printf("torus point: (%0.2f, %0.2f, %0.2f)\n", x, y, z);
			add_point(points, x, y, z);
		}
	}
	return points;
}

/*======== void add_circle() ==========
   Inputs:   struct matrix * edges
            double cx
            double cy
            double r
            double step

   Adds the circle at (cx, cy) with radius r to edges
   ====================*/
void add_circle( struct matrix *edges,
                 double cx, double cy, double cz,
                 double r, int step ) {
	double x0, y0, x1, y1, t;
	int i;

	x0 = r + cx;
	y0 = cy;
	for (i=1; i<=step; i++) {
		t = (double)i/step;
		x1 = r * cos(2 * M_PI * t) + cx;
		y1 = r * sin(2 * M_PI * t) + cy;

		add_edge(edges, x0, y0, cz, x1, y1, cz);
		x0 = x1;
		y0 = y1;
	}
}


/*======== void add_curve() ==========
   Inputs:   struct matrix *edges
         double x0
         double y0
         double x1
         double y1
         double x2
         double y2
         double x3
         double y3
         double step
         int type

   Adds the curve bounded by the 4 points passsed as parameters
   of type specified in type (see matrix.h for curve type constants)
   to the matrix edges
   ====================*/
void add_curve( struct matrix *edges,
                double x0, double y0,
                double x1, double y1,
                double x2, double y2,
                double x3, double y3,
                int step, int type ) {
	double t, x, y;
	int i;
	struct matrix *xcoefs;
	struct matrix *ycoefs;

	xcoefs = generate_curve_coefs(x0, x1, x2, x3, type);
	ycoefs = generate_curve_coefs(y0, y1, y2, y3, type);

	/* print_matrix(xcoefs); */
	/* printf("\n"); */
	/* print_matrix(ycoefs); */

	for (i=1; i<=step; i++) {
		t = (double)i/step;

		x = xcoefs->m[0][0] *t*t*t + xcoefs->m[1][0] *t*t+
		    xcoefs->m[2][0] *t + xcoefs->m[3][0];
		y = ycoefs->m[0][0] *t*t*t + ycoefs->m[1][0] *t*t+
		    ycoefs->m[2][0] *t + ycoefs->m[3][0];

		add_edge(edges, x0, y0, 0, x, y, 0);
		x0 = x;
		y0 = y;
	}

	free_matrix(xcoefs);
	free_matrix(ycoefs);
}



/*======== void add_point() ==========
   Inputs:   struct matrix * points
         int x
         int y
         int z
   Returns:
   adds point (x, y, z) to points and increment points.lastcol
   if points is full, should call grow on points
   ====================*/
void add_point( struct matrix * points, double x, double y, double z) {

	if ( points->lastcol == points->cols )
		grow_matrix( points, points->lastcol + 100 );

	points->m[0][ points->lastcol ] = x;
	points->m[1][ points->lastcol ] = y;
	points->m[2][ points->lastcol ] = z;
	points->m[3][ points->lastcol ] = 1;
	points->lastcol++;
} //end add_point

/*======== void add_edge() ==========
   Inputs:   struct matrix * points
          int x0, int y0, int z0, int x1, int y1, int z1
   Returns:
   add the line connecting (x0, y0, z0) to (x1, y1, z1) to points
   should use add_point
   ====================*/
void add_edge( struct matrix * points,
               double x0, double y0, double z0,
               double x1, double y1, double z1) {
	add_point( points, x0, y0, z0 );
	add_point( points, x1, y1, z1 );
}

/*======== void draw_lines() ==========
   Inputs:   struct matrix * points
         screen s
         color c
   Returns:
   Go through points 2 at a time and call draw_line to add that line
   to the screen
   ====================*/
void draw_lines( struct matrix * points, screen s, zbuffer zb, color c) {

	if ( points->lastcol < 2 ) {
		printf("Need at least 2 points to draw a line!\n");
		return;
	}

	int point;
	for (point=0; point < points->lastcol-1; point+=2)
		draw_line( points->m[0][point],
		           points->m[1][point],
		           points->m[2][point],
		           points->m[0][point+1],
		           points->m[1][point+1],
		           points->m[2][point+1],
		           s, zb, c);
}// end draw_lines


void draw_line(int x0, int y0, double z0,
               int x1, int y1, double z1,
               screen s, zbuffer zb, color c) {



	int xi, yi, xf, yf, xt, yt;
    int d, dx, dy;
    double zi, zt, zf, Dz;


	if(x0 == x1) {
        if(y1 == y0) return;
		if(y0 < y1) {
			yi = y0; yf = y1;
            zi = z0; zf = z1;
		}
		else{
			yi = y1; yf = y0;
            zi = z1; zf = z0;
		}
		xt = x0; yt = yi; zt = zi;
        Dz = (zf - zi) / (yf - yi + 1);

		while(yt <= yf) {
			plot(s, zb, c, xt, yt, zt);
			yt++;
            zt += Dz;
		}
		return;
	}


	else if(x0 < x1) {
		xi = x0; yi = y0;
		xf = x1; yf = y1;

        zi = z0; zf = z1;
	}
	else{
		xi = x1; yi = y1;
		xf = x0; yf = y0;

        zi = z1; zf = z0;
	}
	xt = xi, yt = yi; zt = zi;
    dx = xf - xi, dy = yf - yi;

	int A = dy, B = -dx;
	float m = (double) dy / (double) dx;

	if(m > 0 && m < 1) {
	    //printf("(0,1): %f\n", m);
		d = 2*A + B;
        Dz = (zf - zi) / (xf - xi + 1);

		while(xt <= xf) {
			plot(s, zb, c, xt, yt, zt);
			if(d > 0) {
				d += 2*B;
				yt++;
			}
			d += 2*A;
			xt++;
            zt += Dz;
		}
	}
	else if(m > 1) {
		//printf("(1, +inf): %d %d\n", yf, yi);
		d = A + 2*B;
        Dz = (zf - zi) / (yf - yi + 1);

		while(yt <= yf) {
			plot(s, zb, c, xt, yt, zt);
			if(d < 0) {
				d += 2*A;
				xt++;
			}
			d += 2*B;
			yt++;
            zt += Dz;
		}
	}
	else if(m < 0 && m > -1) {
		//printf("(-1, 0): %f\n", m);
		d = 2*A - B;
        Dz = (zf - zi) / (xf - xi + 1);

		while(xt <= xf) {
			plot(s, zb, c, xt, yt, zt);
			if(d < 0) {
				d -= 2*B;
				yt--;
			}
			d += 2*A;
			xt++;
            zt += Dz;
		}

	}
	else if(m < -1) {
		//printf("(-inf, -1): %f\n", m);
		d = A - 2*B;
        Dz = (zf - zi) / (yf - yi + 1);

		while(yt <= yf) {
			plot(s, zb, c, xt, yt, zt);
			if(d > 0) {
				d += 2*A;
				xt++;
			}
			d -= 2*B;
			yt--;
            zt += Dz;
		}
	}
	else if(m == 0) {
		//printf("{0}: %f\n", m);
        Dz = (zf - zi) / (xf - xi + 1);

		while(xt <= xf) {
			plot(s, zb, c, xt, yt, zt);
			xt++;
            zt += Dz;
		}
	}
	else if(m == 1) {
		//printf("{1}: %f\n", m);
        Dz = (zf - zi) / (xf - xi + 1);

		while(xt <= xf) {
			plot(s, zb, c, xt, yt, zt);
			xt++;
			yt++;
            zt += Dz;
		}
	}
	else if(m == -1) {
	    //printf("{-1}: %f\n", m);
        Dz = (zf - zi) / (xf - xi + 1);

		while(xt <= xf) {
			plot(s, zb, c, xt, yt, zt);
			xt++;
			yt--;
            zt += Dz;
		}

	}
} //end draw_line

package edu.uta.futureye.voronoi;

import java.util.ArrayList;
import java.util.Vector;

/**
 * Voronoi Tessellation by Brute Force Method
 * 
 * @author yuemingl
 *
 */
public class VT2DBruteForce {
	public static double eps = 1e-6;
	public static class Vertex {
		double x;
		double y;
		public Vertex(double x, double y) { this.x = x; this.y = y; }
	}
	
	//A segment starting from a site with angle and length
	public static class Segment {
		double angle;
		double length;
	}
	
	public static class Site {
		double x0;
		double y0;
		ArrayList<Segment> segs = new ArrayList<Segment>(); //sorted by Segment.angle
		boolean isComplete = false; //if the vertices form a polygon, set it to true.
		                            //Indicating that the segments form a 'circle'.
		
		public Vertex getVertex(Segment s1, Segment s2) {
			/**
			 * site: (x0,y0)
			 * 
			 *           line2
			 *            |
			 *  (a,b)     |vertex
			 * --*--------*------------line1
			 *   |        |
			 *  seg1      |
			 *   |        |
			 *   o__seg2__*(c,d)
			 *  (x0,y0)   |
			 *  
			 */
			double a,b,c,d;
			a = x0 + s1.length*Math.cos(s1.angle);
			b = y0 + s1.length*Math.sin(s1.angle);
			c = x0 + s2.length*Math.cos(s2.angle);
			d = y0 + s2.length*Math.sin(s2.angle);
			
			//Vertex = line1 intersects with line2
			if(Math.abs(s1.angle-Math.PI/2) < eps || Math.abs(s1.angle-3*Math.PI/2) < eps) { //s1 is vertical
				//line1 is a horizontal line: y=b 
				if(Math.abs(s2.angle-Math.PI/2) < eps || Math.abs(s2.angle+3*Math.PI/2) < eps) { //s2 is vertical
					return null; //No intersection
				} else if(Math.abs(s2.angle) < eps || Math.abs(s2.angle-Math.PI) < eps){ //s2 is horizontal
					//line2 is a vertical line: x=c
					return new Vertex(c,b);
				} else {
					//line2: y-d=-((c-x0)/(d-y0))*(x-c)
					return new Vertex((b-d)*((d-y0)/(x0-c))+c,b);
				}
			}
			if(Math.abs(s1.angle) < eps || Math.abs(s1.angle-Math.PI) < eps) { //s1 is horizontal
				//line1 is a vertical line: x=a 
				if(Math.abs(s2.angle-Math.PI/2) < eps || Math.abs(s2.angle+3*Math.PI/2) < eps) { //s2 is vertical
					//line2 is a horizontal line: y=d
					return new Vertex(a,d);
				} else if(Math.abs(s2.angle) < eps || Math.abs(s2.angle-Math.PI) < eps){ //s2 is horizontal
					return null; //No intersection
				} else {
					//line2: y-d=-((c-x0)/(d-y0))*(x-c)
					return new Vertex(a, ((x0-c)/(d-y0))*(a-c));
				}
			}
			//General case
			//line1: y-b = -((a-x0)/(b-y0))*(x-a) = k1*(x-a) => y=k1*x+A
			//line2: y-d = -((c-x0)/(d-y0))*(x-c) = k2*(x-c) => y=k2*x+B
			double k1 = (x0-a)/(b-y0);
			double k2 = (x0-c)/(d-y0);
			double A = b-k1*a;
			double B = d-k2*c;
			//k1*x+A==k2*x+B  => x=(B-A)/(k1-k2)
			double x = (B-A)/(k1-k2);
			return new Vertex(x , k1*x+A);
		}

		public void addSegment(Site another) {
			Segment seg = createSegment(another);
			if(segs.size() == 0) { 
				segs.add(seg); 
				return; 
			}
			if(segs.size() == 1) {
				double a = Math.abs(segs.get(0).angle-seg.angle);
				if(a >= eps) {
					if(a <= Math.PI) {
						if(segs.get(0).angle > seg.angle) segs.add(0, seg);
						else segs.add(seg);
					} else {
						if(segs.get(0).angle > seg.angle) segs.add(seg);
						else segs.add(0, seg);
					}
				}
				return;
			}
			
			for(int i=0; i<segs.size()-1; ++i) {
				Segment t = segs.get(i);
				Segment t1 = segs.get(i+1);
				if(Math.abs(seg.angle-t.angle) < eps) {
					//Keep the shorter one
					//TODO: the shorter one may cause other segments be deleted
					//We can sort the neighbor sites by the distance with this site, so this branch will not happen
					if(seg.length < t.length) t.length = seg.length;
				} else if(seg.angle > t.angle && seg.angle < t1.angle) { //这里有问题，例如t.angle=0.01， t2.angle.6.27
					Vertex v = getVertex(t, t1);
					if(v != null) {
						double xx = x0 + seg.length*Math.cos(seg.angle);
						double yy = y0 + seg.length*Math.sin(seg.angle);
						double angle = Math.PI - getAngle(xx-x0,yy-y0, v.x-xx, v.y-yy);
						if(angle > Math.PI/2) //This case, keep the segment, otherwise discard it
							segs.add(i+1,seg);
					} else { //Two segments are on the same line, but have opposite directions
						segs.add(i+1, seg);
					}
					return;
				}
			}
			if(Math.abs(seg.angle-segs.get(segs.size()-1).angle) < eps) {
				//TODO
			}
			Segment head = segs.get(0);
			Segment tail = segs.get(segs.size()-1);
			if(head.angle < tail.angle) {
				
			}
			if(isComplete) {
				if(Math.abs(seg.angle-head.angle) > Math.abs(seg.angle-tail.angle))
					segs.add(0,seg);
				else
					segs.add(seg);
			} else {
				
			}
		}

		private Segment createSegment(Site another) {
			Segment seg = new Segment();
			seg.length = getDistance(this, another);
			if(Math.abs(this.x0-another.x0)<=eps) {
				if(this.y0 < another.y0)
					seg.angle = Math.PI/2;
				else 
					seg.angle = -Math.PI/2;
			} else {
				seg.angle = Math.atan2(another.y0-this.y0, another.x0-this.x0);
			}
			return seg;
		}
		
		public static double getDistance(Site s1, Site s2) {
			return Math.sqrt((s1.x0-s2.x0)*(s1.x0-s2.x0)+(s1.y0-s2.y0)*(s1.y0-s2.y0));
		}
		
		//Get the angle between vector u and vector v
		public static double getAngle(double ux, double uy, double vx, double vy) {
			//u.v=|u|*|v|*cos(angle)
			double lu = Math.sqrt(ux*ux+uy*uy);
			double lv = Math.sqrt(vx*vx+vy*vy);
			double dot = ux*vx+uy*vy;
			return Math.acos(dot/(lu*lv));
		}
		
	}
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub

	}

}

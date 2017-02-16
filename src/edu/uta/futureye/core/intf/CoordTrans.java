package edu.uta.futureye.core.intf;

import java.util.HashMap;

import edu.uta.futureye.function.intf.MathFunc;

/**
 * Coordinate transformation
 * 
 */
public interface CoordTrans {
	/**
	 * Return the coordinates.
	 * For example, the area coordinates on a triangle is r,s,t
	 * where t=1-r-s
	 * @return
	 */
	public MathFunc[] getCoords();
	
	/**
	 * Return the coordinate transform map
	 * For example,
	 *  x = x1*r + x2*s + x3*(1-r-s);
	 *  y = y1*r + y2*s + y3*(1-r-s);
	 *  HashMap<String, MathFunc> map = new HashMap<>();
	 *  map.put("x", x);
	 *  map.put("y", y);
	 * @return
	 */
	public HashMap<String, MathFunc> getCoordTransMap();

	/**
	 * Return the Jacobian of the coordinate transformation
	 * For example, 2D coordinate transform
	 * det(x_r, x_s)
	 *    (y_r, y_s)
	 * @return
	 */
	public MathFunc getJacobian();
}

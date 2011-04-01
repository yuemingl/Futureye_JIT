package edu.uta.futureye.core.geometry;

/**
 * 有限元空间点
 * 
 * @author liuyueming
 *
 */
public interface Point extends GeoEntity {
	
	int dim();
	
	double coord(int index);
	
	double[] coords();
	
	void setCoord(int index,double value);
	
	/**
	 * 判断坐标点是否相等
	 * @param p
	 * @return
	 */
	boolean coordEquals(Point p);
	
	/**
	 * 获取索引编号，
	 * 对于全局结点，返回全局编号，
	 * 对于局部结点，返回局部编号
	 * 
	 * @return
	 */
	int getIndex();
}

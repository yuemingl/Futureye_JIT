package edu.uta.futureye.core.intf;

public interface Point extends GeoEntity {
	double coord(int index);
	int dim();
	public boolean coordEquals(Point p);
}

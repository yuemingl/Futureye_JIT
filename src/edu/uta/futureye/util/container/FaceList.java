package edu.uta.futureye.util.container;

import edu.uta.futureye.core.Face;

/**
 * Face List Class
 * （全局）面类列表
 * 
 * @author liuyueming
 *
 */
public class FaceList extends ObjList<Face> {
	@Override
	public String toString() {
		return "FaceList"+objs.toString();
	}
}

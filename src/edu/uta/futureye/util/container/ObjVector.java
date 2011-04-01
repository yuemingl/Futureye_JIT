package edu.uta.futureye.util.container;

import java.util.List;
import java.util.Vector;

import edu.uta.futureye.util.FutureyeException;

/**
 * Object vector container
 * * index starts from 1
 * * <tt>null</tt> element is allowed.
 * * size can be specified.
 * 
 * @author liuyueming
 *
 * @param <T>
 */
public class ObjVector<T> {
	protected Vector<T> objs = new Vector<T>();
	
	public ObjVector() {
	}
	
	public ObjVector(int size) {
		objs.setSize(size);
	}
	
	public ObjVector(T ...es) {
		for(T e : es) this.add(e);
	}

	/**
	 * @param index starts from 1,2,3...
	 * @return
	 */
	public T at(int index) {
		if(index < 1)
			throw new FutureyeException("ERROR: ObjVector index should be >=1, index="+index);
		if(objs.size() < index)
			return null;
		return objs.get(index-1);
	}

	public ObjVector<T> add(T e) {
		this.objs.add(e);
		return this;
	}
	
	public ObjVector<T> addAll(ObjList<T> list) {
		if(list == null) return this;
		this.objs.addAll(list.objs);
		return this;
	}
	
	/**
     * Replaces the element at the specified position in this ObjVector with 
     * the specified element <tt>e</tt>.
     * 
	 * @param index
	 * @param e
	 * @return the whole new ObjVector
	 */
	public ObjVector<T> set(int index,T e) {
		if(index < 1)
			throw new FutureyeException("ERROR: ObjVector index should be >=1, index="+index);
		this.objs.set(index-1, e);
		return this;
	}
	
	public int size() {
		return objs.size();
	}
	
    /**
     * Sets the size of this ObjVector. 
     * If the new size is greater than the current size, 
     * new <tt>null</tt> items are added to the end of the vector. 
     * If the new size is less than the current size, 
     * all components at index <tt>newSize+1</tt> and greater are discarded.
     *
     * @param  newSize   the new size of this vector
     */
	public void setSize(int newSize) {
		objs.setSize(newSize);
	}
	
	public void clear() {
		objs.clear();
	}
	
	public T remove(int index) {
		return objs.remove(index-1);
	}
	
	public boolean remove(T e) {
		return objs.remove(e);
	}
	
	public ObjVector<T> subVector(int begin,int end) {
		ObjVector<T> rlt = new ObjVector<T>();
		for(int i=begin;i<=end;i++)
			rlt.add(this.at(i));
		return rlt;
	}
	
	public ObjVector<T> subVector(int begin,int end,int step) {
		ObjVector<T> rlt = new ObjVector<T>();
		for(int i=begin;i<=end;i+=step)
			rlt.add(this.at(i));
		return rlt;
	}
	
	public ObjVector<T> subVector(List<Integer> set) {
		ObjVector<T> rlt = new ObjVector<T>();
		for(Integer i : set)
			rlt.add(this.at(i));
		return rlt;
	}
	
	public List<T> toList() {
		return objs;
	}
	public Object[] toArray() {
		return objs.toArray();
	}
	
	/**
	 * 从一个下标为0,1,2,...的list<T>构建
	 * 为一个下标为1,2,3,...的ObjList<T>
	 * @param list
	 */
	public ObjVector<T> fromList(List<T> list) {
		objs.clear();
		objs.addAll(list);
		return this;
	}
	public ObjVector<T> fromArray(T[] array) {
		objs.clear();
		for(int i=0;i<array.length;i++)
			objs.add((T)array[i]);
		return this;
	}
	
	public String toString() {
		return objs.toString();
	}
}

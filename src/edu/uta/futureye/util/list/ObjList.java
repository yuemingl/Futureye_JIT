package edu.uta.futureye.util.list;

import java.util.ArrayList;
import java.util.List;

import edu.uta.futureye.util.FutureyeException;

public class ObjList<T> {
	protected List<T> objs = new ArrayList<T>();
	
	public ObjList() {
	}
	
	public ObjList(T ...es) {
		for(T e : es) this.add(e);
	}

	/**
	 * @param index start from 1,2,3...
	 * @return
	 */
	public T at(int index) {
		if(index < 1) {
			FutureyeException e = new FutureyeException("ERROR: ObjList index="+index);
			e.printStackTrace();
			System.exit(0);
		}
		return objs.get(index-1);
	}

	public void add(T e) {
		this.objs.add(e);
	}
	
	public void addAll(ObjList<T> list) {
		if(list == null) return;
		this.objs.addAll(list.objs);
	}
	
	public int size() {
		return objs.size();
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
	
	public ObjList<T> subList(int begin,int end) {
		ObjList<T> rlt = new ObjList<T>();
		for(int i=begin;i<=end;i++)
			rlt.add(this.at(i));
		return rlt;
	}
	
	public ObjList<T> subList(List<Integer> set) {
		ObjList<T> rlt = new ObjList<T>();
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
	public ObjList<T> fromList(List<T> list) {
		objs.clear();
		objs.addAll(list);
		return this;
	}
	public ObjList<T> fromArray(T[] array) {
		objs.clear();
		for(int i=0;i<array.length;i++)
			objs.add((T)array[i]);
		return this;
	}
	
	public String toString() {
		return objs.toString();
	}
}

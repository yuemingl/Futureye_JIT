package edu.uta.futureye.util;

public class FutureyeException extends RuntimeException {

	private static final long serialVersionUID = 1L;

	public FutureyeException() {
		super("");
	}

	public FutureyeException(String info) {
		super(info);
	}

}

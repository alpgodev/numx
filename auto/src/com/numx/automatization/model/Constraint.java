package com.numx.automatization.model;

/**
 * DOCUMENT ME!
 */
public class Constraint implements Comparable {
   
	private String _constraintType;
	private String _constraintValue;
	private String _constraintValue2;
	private String _constraintValueJava;
	private String _constraintValueJava2;
	private String _constraintFirstIndex;
	private String _constraintLastIndex;
	private String _constraintFirstIndexJava;
	private String _constraintLastIndexJava;
	

    /**
     * Creates a new Parameter object.
     *
     * @constraint type DOCUMENT ME!
     * @constraint value DOCUMENT ME!
     */
    public Constraint(String type, String value, String value2, String valueJava, 
    		String valueJava2, String firstIndex, String lastIndex, 
    		String firstIndexJava, String lastIndexJava) {
        _constraintType = type;
        _constraintValue = value;
		_constraintValue2 = value2;
		_constraintValueJava = valueJava;
		_constraintValueJava2 = valueJava2;
		_constraintFirstIndex = firstIndex;
		_constraintLastIndex = lastIndex;
		_constraintFirstIndexJava = firstIndexJava;
		_constraintLastIndexJava = lastIndexJava;
    }

    /**
     * DOCUMENT ME!
     *
     * @return DOCUMENT ME!
     */
    public String get_constraintType() {
        return _constraintType;
    }

    /**
     * DOCUMENT ME!
     *
     * @constraint columns DOCUMENT ME!
     */
    public void set_constraintType(String type) {
        _constraintType = type;
    }

    /**
     * DOCUMENT ME!
     *
     * @return DOCUMENT ME!
     */
    public String get_firstIndexJava() {
        return _constraintFirstIndexJava;
    }

    /**
     * DOCUMENT ME!
     *
     * @constraint value DOCUMENT ME!
     */
    public void set_firstIndexJava(String firstIndexJava) {
    	_constraintFirstIndexJava = firstIndexJava;
    }


    /**
     * DOCUMENT ME!
     *
     * @return DOCUMENT ME!
     */
    public String get_lastIndexJava() {
        return _constraintLastIndexJava;
    }

    /**
     * DOCUMENT ME!
     *
     * @constraint value DOCUMENT ME!
     */
    public void set_lastIndexJava(String lastIndexJava) {
    	_constraintLastIndexJava = lastIndexJava;
    }


    /**
     * DOCUMENT ME!
     *
     * @return DOCUMENT ME!
     */
    public String get_firstIndex() {
        return _constraintFirstIndex;
    }

    /**
     * DOCUMENT ME!
     *
     * @constraint value DOCUMENT ME!
     */
    public void set_firstIndex(String firstIndex) {
    	_constraintFirstIndex = firstIndex;
    }


    /**
     * DOCUMENT ME!
     *
     * @return DOCUMENT ME!
     */
    public String get_lastIndex() {
        return _constraintLastIndex;
    }

    /**
     * DOCUMENT ME!
     *
     * @constraint value DOCUMENT ME!
     */
    public void set_lastIndex(String lastIndex) {
    	_constraintLastIndex = lastIndex;
    }

    /**
     * DOCUMENT ME!
     *
     * @return DOCUMENT ME!
     */
    public String get_constraintValue() {
        return _constraintValue;
    }

    /**
     * DOCUMENT ME!
     *
     * @constraint value DOCUMENT ME!
     */
    public void set_constraintValue(String value) {
        _constraintValue = value;
    }

	/**
	 * DOCUMENT ME!
	 *
	 * @return DOCUMENT ME!
	 */
	public String get_constraintValue2() 
	{
		return _constraintValue2;
	}

	/**
	 * DOCUMENT ME!
	 *
	 * @constraint value DOCUMENT ME!
	 */
	public void set_constraintValue2(String value2) 
	{
		_constraintValue2 = value2;
	}
	
	/**
	 * DOCUMENT ME!
	 *
	 * @return DOCUMENT ME!
	 */
	public String get_constraintValueJava() 
	{
		return _constraintValueJava;
	}

	/**
	 * DOCUMENT ME!
	 *
	 * @constraint value DOCUMENT ME!
	 */
	public void set_constraintValueJava(String valueJava) 
	{
		_constraintValueJava = valueJava;
	}

	/**
	 * DOCUMENT ME!
	 *
	 * @return DOCUMENT ME!
	 */
	public String get_constraintValueJava2() 
	{
		return _constraintValueJava2;
	}

	/**
	 * DOCUMENT ME!
	 *
	 * @constraint value DOCUMENT ME!
	 */
	public void set_constraintValueJava2(String valueJava2) 
	{
		_constraintValueJava2 = valueJava2;
	}

	public int compareTo(Object arg0) 
	{
		return 0;//_module.compareTo(((Constraint) arg0).get_moduleName());	
	}


}

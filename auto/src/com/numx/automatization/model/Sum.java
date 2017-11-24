package com.numx.automatization.model;

/**
 * DOCUMENT ME!
 */
public class Sum implements Comparable {
   
	private String _sumName;
	private String _sumType;
	private String _sumValue;
	private String _condValue1;
	private String _condValue2;
	private String _sumValueJNI;
	private String _condValue1JNI;
	private String _condValue2JNI;
	private String _sumValueVB;
	private String _condValue1VB;
	private String _condValue2VB;
	private String _condType;

    /**
     * Creates a new Parameter object.
     *
     * @sum type DOCUMENT ME!
     * @sum value DOCUMENT ME!
     */
	public Sum(String name, String type, String value, String condvalue1, String condvalue2, String condtype,
		String valueJNI, String condvalue1JNI, String condvalue2JNI, String valueVB, String condValue1VB,
		String condValue2VB) {
		_sumName = name;
        _sumType = type; // normal, equalCond
        _sumValue = value;
		_condValue1 = condvalue1;
		_condValue2 = condvalue2;
		_sumValueJNI = valueJNI;
		_condValue1JNI = condvalue1JNI;
		_condValue2JNI = condvalue2JNI;
		_sumValueVB = valueVB;
		_condValue1VB = condValue1VB;
		_condValue2VB = condValue2VB;
		_condType = condtype;
    }
	
	/**
	 * DOCUMENT ME!
	 *
	 * @return DOCUMENT ME!
	 */
	public String get_sumName() 
	{
		return _sumName;
	}

	/**
	 * DOCUMENT ME!
	 *
	 * @sum columns DOCUMENT ME!
	 */
	public void set_sumName(String name) 
	{
		_sumName = name;
	}

	
	/**
	 * DOCUMENT ME!
	 *
	 * @return DOCUMENT ME!
	 */
	public String get_condValue1() 
	{
		return _condValue1;
	}

	/**
	 * DOCUMENT ME!
	 *
	 * @sum columns DOCUMENT ME!
	 */
	public void set_condValue1(String condvalue1) 
	{
		_condValue1 = condvalue1;
	}


	/**
	 * DOCUMENT ME!
	 *
	 * @return DOCUMENT ME!
	 */
	public String get_condValue2() 
	{
		return _condValue2;
	}

	/**
	 * DOCUMENT ME!
	 *
	 * @sum columns DOCUMENT ME!
	 */
	public void set_condValue2(String condvalue2) 
	{
		_condValue2 = condvalue2;
	}
	
	/**
	 * DOCUMENT ME!
	 *
	 * @sum columns DOCUMENT ME!
	 */
	public String get_condType() 
	{
		return _condType;
	}



	/**
	 * DOCUMENT ME!
	 *
	 * @sum columns DOCUMENT ME!
	 */
	public void set_condType(String condtype) 
	{
		_condType = condtype;
	}

    /**
     * DOCUMENT ME!
     *
     * @return DOCUMENT ME!
     */
    public String get_sumType() 
	{
        return _sumType;
    }

    /**
     * DOCUMENT ME!
     *
     * @sum columns DOCUMENT ME!
     */
    public void set_sumType(String type) {
        _sumType = type;
    }

    /**
     * DOCUMENT ME!
     *
     * @return DOCUMENT ME!
     */
    public String get_sumValue() {
        return _sumValue;
    }

    /**
     * DOCUMENT ME!
     *
     * @sum value DOCUMENT ME!
     */
    public void set_sumValue(String value) {
        _sumValue = value;
    }

	/**
	 * DOCUMENT ME!
	 *
	 * @return DOCUMENT ME!
	 */
	public String get_condValue1JNI() 
	{
		return _condValue1JNI;
	}

	/**
	 * DOCUMENT ME!
	 *
	 * @sum columns DOCUMENT ME!
	 */
	public void set_condValue1JNI(String condvalue1JNI) 
	{
		_condValue1JNI = condvalue1JNI;
	}


	/**
	 * DOCUMENT ME!
	 *
	 * @return DOCUMENT ME!
	 */
	public String get_condValue2JNI() 
	{
		return _condValue2JNI;
	}

	/**
	 * DOCUMENT ME!
	 *
	 * @sum columns DOCUMENT ME!
	 */
	public void set_condValue2JNI(String condvalue2JNI) 
	{
		_condValue2JNI = condvalue2JNI;
	}
	
	/**
	 * DOCUMENT ME!
	 *
	 * @return DOCUMENT ME!
	 */
	public String get_sumValueJNI() 
	{
		return _sumValueJNI;
	}

	/**
	 * DOCUMENT ME!
	 *
	 * @sum value DOCUMENT ME!
	 */
	public void set_sumValueJNI(String valueJNI) 
	{
		_sumValueJNI = valueJNI;
	}

	
	/**
	 * DOCUMENT ME!
	 *
	 * @return DOCUMENT ME!
	 */
	public String get_condValue1VB() 
	{
		return _condValue1VB;
	}

	/**
	 * DOCUMENT ME!
	 *
	 * @sum columns DOCUMENT ME!
	 */
	public void set_condValue1VB(String condvalue1VB) 
	{
		_condValue1VB = condvalue1VB;
	}


	/**
	 * DOCUMENT ME!
	 *
	 * @return DOCUMENT ME!
	 */
	public String get_condValue2VB() 
	{
		return _condValue2VB;
	}

	/**
	 * DOCUMENT ME!
	 *
	 * @sum columns DOCUMENT ME!
	 */
	public void set_condValue2VB(String condvalue2VB) 
	{
		_condValue2VB = condvalue2VB;
	}
	
	/**
	 * DOCUMENT ME!
	 *
	 * @return DOCUMENT ME!
	 */
	public String get_sumValueVB() 
	{
		return _sumValueVB;
	}

	/**
	 * DOCUMENT ME!
	 *
	 * @sum value DOCUMENT ME!
	 */
	public void set_sumValueVB(String valueVB) 
	{
		_sumValueVB = valueVB;
	}



	/**
	 * DOCUMENT ME!
	 *
	 * @sum columns DOCUMENT ME!
	 */
	public int compareTo(Object arg0) 
	{
		return 0;//_module.compareTo(((sum) arg0).get_moduleName());	
	}


}

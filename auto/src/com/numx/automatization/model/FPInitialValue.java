package com.numx.automatization.model;


/**
 * DOCUMENT ME!
 */
public class FPInitialValue implements Comparable 
{
   
	private String _fpInitialValueValue;
	private String _fpInitialValueIndex;

    /**
     * Creates a new Parameter object.
     *
     * @sum name DOCUMENT ME!
     * @sum parameters DOCUMENT ME!
     */
	public FPInitialValue(String value, String index) {
		_fpInitialValueValue = value;
		_fpInitialValueIndex = index;
    }
	
	/**
	 * DOCUMENT ME!
	 *
	 * @return DOCUMENT ME!
	 */
	public String get_fpInitialValueValue(){
		return _fpInitialValueValue;
	}

	/**
	 * DOCUMENT ME!
	 *
	 * @sum columns DOCUMENT ME!
	 */
	public void set_fpInitialValueValue(String value){
		_fpInitialValueValue = value;
	}

	/**
	 * DOCUMENT ME!
	 *
	 * @return DOCUMENT ME!
	 */
	public String get_fpInitialValueIndex()
	{
		return _fpInitialValueIndex;
	}

	/**
	 * DOCUMENT ME!
	 *
	 * @sum columns DOCUMENT ME!
	 */
	public void set_fpInitialValueIndex(String index)
	{
		_fpInitialValueIndex = index;
	}

	public int compareTo(Object arg0) 
	{
		return 0;//_module.compareTo(((sum) arg0).get_moduleName());	
	}


}

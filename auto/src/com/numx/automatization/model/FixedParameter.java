package com.numx.automatization.model;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * DOCUMENT ME!
 */
public class FixedParameter implements Comparable 
{
   
	private String _fixedParameterName;
	private String _fixedParameterValue;
	private String _fixedParameterType;
	private List fpInitialValues;

    /**
     * Creates a new Parameter object.
     *
     * @sum name DOCUMENT ME!
     * @sum parameters DOCUMENT ME!
     */
	public FixedParameter(String name, String value, String type) {
		_fixedParameterName = name;
		_fixedParameterValue = value;
		_fixedParameterType = type;
		fpInitialValues = new ArrayList();
    }
	
	/**
	 * DOCUMENT ME!
	 *
	 * @return DOCUMENT ME!
	 */
	public String get_fixedParameterName(){
		return _fixedParameterName;
	}

	/**
	 * DOCUMENT ME!
	 *
	 * @sum columns DOCUMENT ME!
	 */
	public void set_fixedParameterName(String name){
		_fixedParameterName = name;
	}

	/**
	 * DOCUMENT ME!
	 *
	 * @return DOCUMENT ME!
	 */
	public String get_fixedParameterValue()
	{
		return _fixedParameterValue;
	}

	/**
	 * DOCUMENT ME!
	 *
	 * @sum columns DOCUMENT ME!
	 */
	public void set_fixedParameterValue(String value)
	{
		_fixedParameterValue = value;
	}

	/**
	 * DOCUMENT ME!
	 *
	 * @return DOCUMENT ME!
	 */
	public String get_fixedParameterType()
	{
		return _fixedParameterType;
	}

	/**
	 * DOCUMENT ME!
	 *
	 * @sum columns DOCUMENT ME!
	 */
	public void set_fixedParameterType(String type)
	{
		_fixedParameterType = type;
	}

	/**
	 * Get parameters ordered by position
	 * @return
	 */
	public List getFPInitialValues() 
	{
		return fpInitialValues;
	}

	public void setFPInitialValues(List fpinitialvalues) 
	{
		this.fpInitialValues = fpinitialvalues;
	}

	public void addFPInitialValue(FPInitialValue fpinitialvalue) 
	{
		fpInitialValues.add(fpinitialvalue);
		Collections.sort(fpInitialValues);
	}


	public int compareTo(Object arg0) 
	{
		return 0;//_module.compareTo(((sum) arg0).get_moduleName());	
	}


}

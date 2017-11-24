package com.numx.automatization.model;


/**
 * DOCUMENT ME!
 */
public class GetWorkspace implements Comparable 
{
   
	private String _gwName;
	private String _gwParamType;
	private String _gwType;// all double, int or no.

    /**
     * Creates a new Parameter object.
     *
     * @sum name DOCUMENT ME!
     * @sum parameters DOCUMENT ME!
     */
	public GetWorkspace(String name, String type, String paramType) {
		_gwName = name;
		_gwType = type;
		_gwParamType = paramType;
    }
	
	/**
	 * DOCUMENT ME!
	 *
	 * @return DOCUMENT ME!
	 */
	public String get_gwName()
	{
		return _gwName;
	}

	/**
	 * DOCUMENT ME!
	 *
	 * @sum columns DOCUMENT ME!
	 */
	public void set_gwName(String name){
		_gwName = name;
	}
	
	/**
	 * DOCUMENT ME!
	 *
	 * @return DOCUMENT ME!
	 */
	public String get_gwParamType()
	{
		return _gwParamType;
	}

	/**
	 * DOCUMENT ME!
	 *
	 * @sum columns DOCUMENT ME!
	 */
	public void set_gwParamType(String paramType)
	{
		_gwParamType = paramType;
	}

	/**
	 * DOCUMENT ME!
	 *
	 * @return DOCUMENT ME!
	 */
	public String get_gwType()
	{
		return _gwType;
	}

	/**
	 * DOCUMENT ME!
	 *
	 * @sum columns DOCUMENT ME!
	 */
	public void set_gwType(String type)
	{
		_gwType = type;
	}



	public int compareTo(Object arg0) 
	{
		return 0;//_module.compareTo(((sum) arg0).get_moduleName());	
	}


}

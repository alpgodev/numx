package com.numx.automatization.model;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * DOCUMENT ME!
 */
public class Parameter implements Comparable 
{
    /**
     * DOCUMENT ME!
     */
    public final static String SIMEXT = "void (*F77_FUNCTION(simext))(int*, int*, double*, int*, double*, int*, int*, int*, double*, double*, double*, int*)";

    /**
     * DOCUMENT ME!
     */
    public final static String SIMUL = "void (*F77_FUNCTION(simul))(int*, void (*F77_FUNCTION(simext))(int*, int*, double*, int*, double*, int*, int*, int*, double*, double*, double*, int*), int*, double*, double*, double*, int*, double*)";

    /**
     * DOCUMENT ME!
     */
    public final static String USERPO = "void (*F77_FUNCTION(userpo))(int*, int*, double*, double*, int*, double*, int*, double*, double*, int*)";
   
    private int _parameterId;
    private String _parameterName;
    private String _parameterType;
    private String _parameterDescription = "";
    private boolean _inputParameter;
    private boolean _lastInput;
    private boolean _outFunction;
    private int _parameterOrder;
    private String _nbRows;
	private String _nbColumns;
	private String _nbTestRows;
	private String _nbTestColumns;
    private boolean _zeroLength;
    private boolean _random;
	private List constraints;
	private List sums;

    /**
     * Creates a new Parameter object.
     *
     * @param id DOCUMENT ME!
     * @param name DOCUMENT ME!
     * @param type DOCUMENT ME!
     * @param description DOCUMENT ME!
     * @param parameter DOCUMENT ME!
     * @param order DOCUMENT ME!
     * @param rows DOCUMENT ME!
     * @param columns DOCUMENT ME!
     * @param length DOCUMENT ME!
     * @param _random DOCUMENT ME!
     */
    public Parameter(int id, String name, String type, String description,
        boolean parameter, boolean lastInput, boolean outfunction, int order, String rows, String columns, String testRows, String testColumns,
        boolean length, boolean random) {
        _parameterId = id;
        _parameterName = name;
        _parameterType = type;
        _parameterDescription = description;
        _inputParameter = parameter;
        _lastInput = lastInput;
        _outFunction = outfunction;
        _parameterOrder = order;
        _nbRows = rows;
        _nbColumns = columns;
		_nbTestRows = testRows;
		_nbTestColumns = testColumns;
        _zeroLength = length;
        _random = random;
		constraints = new ArrayList();
		sums = new ArrayList();
    }

    /**
     * Creates a new Parameter object.
     *
     * @param id DOCUMENT ME!
     * @param name DOCUMENT ME!
     * @param type DOCUMENT ME!
     * @param description DOCUMENT ME!
     * @param parameter DOCUMENT ME!
     * @param order DOCUMENT ME!
     * @param rows DOCUMENT ME!
     * @param columns DOCUMENT ME!
     */
  /*  public Parameter(int id, String name, String type, String description,
        boolean parameter, int order, String rows, String columns) {
        _parameterId = id;
        _parameterName = name;
        _parameterType = type;
        _parameterDescription = description;
        _inputParameter = parameter;
        _parameterOrder = order;
        _nbRows = rows;
        _nbColumns = columns;
    }*/
    
    /**
     * DOCUMENT ME!
     *
     * @return DOCUMENT ME!
     */
    public boolean is_lastInput() {
        return _lastInput;
    }

    /**
     * DOCUMENT ME!
     *
     * @param parameter DOCUMENT ME!
     */
    public void set_lastInput(boolean lastInput) {
        _lastInput = lastInput;
    }
    
    
    /**
     * DOCUMENT ME!
     *
     * @return DOCUMENT ME!
     */
    public boolean is_inputParameter() {
        return _inputParameter;
    }

    /**
     * DOCUMENT ME!
     *
     * @param parameter DOCUMENT ME!
     */
    public void set_inputParameter(boolean parameter) {
        _inputParameter = parameter;
    }
    
    /**
     * DOCUMENT ME!
     *
     * @return DOCUMENT ME!
     */
    public boolean is_outFunction() {
        return _outFunction;
    }

    /**
     * DOCUMENT ME!
     *
     * @param parameter DOCUMENT ME!
     */
    public void set_outFunction(boolean outfunction) {
    	_outFunction = outfunction;
    }
    
    
    /**
     * DOCUMENT ME!
     *
     * @return DOCUMENT ME!
     */
    public String get_nbColumns() {
        return _nbColumns;
    }

    /**
     * DOCUMENT ME!
     *
     * @param columns DOCUMENT ME!
     */
    public void set_nbColumns(String columns) {
        _nbColumns = columns;
    }

    /**
     * DOCUMENT ME!
     *
     * @return DOCUMENT ME!
     */
    public String get_nbRows() {
        return _nbRows;
    }

    /**
     * DOCUMENT ME!
     *
     * @param rows DOCUMENT ME!
     */
    public void set_nbRows(String rows) {
        _nbRows = rows;
    }
	
	/**
	 * DOCUMENT ME!
	 *
	 * @return DOCUMENT ME!
	 */
	public String get_nbTestColumns() 
	{
		return _nbTestColumns;
	}

	/**
	 * DOCUMENT ME!
	 *
	 * @param columns DOCUMENT ME!
	 */
	public void set_nbTestColumns(String testcolumns) 
	{
		_nbTestColumns = testcolumns;
	}

	/**
	 * DOCUMENT ME!
	 *
	 * @return DOCUMENT ME!
	 */
	public String get_nbTestRows() 
	{
		return _nbTestRows;
	}

	/**
	 * DOCUMENT ME!
	 *
	 * @param rows DOCUMENT ME!
	 */
	public void set_nbTestRows(String testrows) 
	{
		_nbTestRows = testrows;
	}

    /**
     * DOCUMENT ME!
     *
     * @return DOCUMENT ME!
     */
    public String get_parameterDescription() 
	{
        return _parameterDescription;
    }

    /**
     * DOCUMENT ME!
     *
     * @param description DOCUMENT ME!
     */
    public void set_parameterDescription(String description) {
        _parameterDescription = description;
    }

    /**
     * DOCUMENT ME!
     *
     * @return DOCUMENT ME!
     */
    public String get_parameterName() {
        return _parameterName;
    }

    /**
     * DOCUMENT ME!
     *
     * @param name DOCUMENT ME!
     */
    public void set_parameterName(String name) {
        _parameterName = name;
    }

    /**
     * DOCUMENT ME!
     *
     * @return DOCUMENT ME!
     */
    public String getFunctionHeaderByParameterName() {
        if (_parameterName.equals("simext")) {
            return SIMEXT;
        } else if (_parameterName.equals("simul")) {
            return SIMUL;
        } else if (_parameterName.equals("userpo")) {
            return USERPO;
        } else {
            return "";
        }
    }

    /**
     * DOCUMENT ME!
     *
     * @return DOCUMENT ME!
     */
    public int get_parameterOrder() {
        return _parameterOrder;
    }

    /**
     * DOCUMENT ME!
     *
     * @param order DOCUMENT ME!
     */
    public void set_parameterOrder(int order) {
        _parameterOrder = order;
    }

    /**
     * DOCUMENT ME!
     *
     * @return DOCUMENT ME!
     */
    public String get_parameterType() {
        return _parameterType;
    }

    /**
     * DOCUMENT ME!
     *
     * @param type DOCUMENT ME!
     */
    public void set_parameterType(String type) {
        _parameterType = type;
    }

    /**
     * DOCUMENT ME!
     *
     * @return DOCUMENT ME!
     */
    public int get_parameterId() {
        return _parameterId;
    }

    /**
     * DOCUMENT ME!
     *
     * @param id DOCUMENT ME!
     */
    public void set_parameterId(int id) {
        _parameterId = id;
    }

    /**
     * DOCUMENT ME!
     *
     * @param arg0 DOCUMENT ME!
     *
     * @return DOCUMENT ME!
     */
    public int compareTo(Object arg0) {
        return _parameterOrder - ((Parameter) arg0).get_parameterOrder();
    }

    /**
     * DOCUMENT ME!
     *
     * @return DOCUMENT ME!
     */
    public boolean is_random() {
        return _random;
    }

    /**
     * DOCUMENT ME!
     *
     * @param _random DOCUMENT ME!
     */
    public void set_random(boolean _random) {
        this._random = _random;
    }

    /**
     * DOCUMENT ME!
     *
     * @return DOCUMENT ME!
     */
    public boolean is_zeroLength() {
        return _zeroLength;
    }

    /**
     * DOCUMENT ME!
     *
     * @param length DOCUMENT ME!
     */
    public void set_zeroLength(boolean length) {
        _zeroLength = length;
    }


	/**
	 * Get functions ordered by name
	 * @return
	 */
	public List get_Constraints() 
	{
		return constraints;
	}

	public void set_Constraints(List constraints) 
	{
		this.constraints = constraints;
	}

	public void addConstraint(Constraint constraint) 
	{
		constraints.add(constraint);
		Collections.sort(constraints);
	}

	/**
	 * Get functions ordered by name
	 * @return
	 */
	public List get_Sums() 
	{
		return sums;
	}

	public void set_Sums(List sums) 
	{
		this.sums = sums;
	}

	public void addSum(Sum sum) 
	{
		sums.add(sum);
		Collections.sort(sums);
	}
}

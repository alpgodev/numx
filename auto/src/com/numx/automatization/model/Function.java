package com.numx.automatization.model;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;


public class Function implements Comparable {
	
	private int _functionId;
	private String _functionName;
	private String _functionShortDescription = "";
	private String _functionDescription = "";
	private boolean _onlyC;
	private boolean _release;
	private boolean _istest;
	private boolean _isFortran;
	private String _fortranFunction;
	private List getWorkspaces; 
	private List parameters;
	private List fixedParameters;
	private String _gwFunctionName="";
	
	public Function(int id, String name, String shortDescription, String description, boolean onlyc, boolean release, boolean istest, boolean isfortran,
		String function, String gwFunctionName) {
		_functionId = id;
		_functionName = name;
		_functionShortDescription = shortDescription;
		_functionDescription = description;
		_onlyC = onlyc;
		_release = release;
		_istest = istest;
		_fortranFunction = function;
		getWorkspaces = new ArrayList();
		parameters = new ArrayList();
		fixedParameters = new ArrayList();
		_gwFunctionName = gwFunctionName;
		_isFortran= isfortran;
	}

	public String get_functionDescription() {
		return _functionDescription;
	}

	public void set_functionDescription(String description) {
		_functionDescription = description;
	}

	public int get_functionId() {
		return _functionId;
	}

	public void set_functionId(int id) {
		_functionId = id;
	}

	public String get_functionName() {
		return _functionName;
	}

	public void set_functionName(String name) {
		_functionName = name;
	}

	public String get_functionShortDescription() {
		return _functionShortDescription;
	}

	public void set_functionShortDescription(String shortDescription) {
		_functionShortDescription = shortDescription;
	}

	public boolean is_onlyC() {
		return _onlyC;
	}

	public void set_onlyC(boolean _onlyc) {
		_onlyC = _onlyc;
	}	

	public boolean is_fortran() 
	{
		return _isFortran;
	}

	public void set_isFortran(boolean _isfortran) 
	{
		_isFortran = _isfortran;
	}	



	public boolean is_test() 
	{
		return _istest;
	}

	public void set_istest(boolean _isTest) 
	{
		_istest = _isTest;
	}
	
	
	public String get_gwFunctionName() 
	{
		return _gwFunctionName;
	}

	public void set_gwFunctionName(String gwfunctionname) 
	{
		_gwFunctionName = gwfunctionname;
	}

	public boolean gwFunctionName_is_defined() 
	{
		if (_gwFunctionName == "")
		{
			return false;
		}
		else
		{
			return true;
		}
		
	}

	/**
	 * Get parameters ordered by position
	 * @return
	 */
	public List getParameters() 
	{
		return parameters;
	}

	public void setParameters(List parameters) {
		this.parameters = parameters;
	}

	public void addParameter(Parameter parameter) {
		parameters.add(parameter);
		Collections.sort(parameters);
	}

	/**
	 * Get parameters ordered by position
	 * @return
	 */
	public List getGetWorkspaces() 
	{
		return getWorkspaces;
	}

	public void setGetWorkspaces(List getworkspace) 
	{
		this.getWorkspaces = getworkspace;
	}

	public void addGetWorkspace(GetWorkspace getworkspace) 
	{
		getWorkspaces.add(getworkspace);
		Collections.sort(getWorkspaces);
	}
	
	/**
	 * Get parameters ordered by position
	 * @return
	 */
	public List getFixedParameter() 
	{
		return fixedParameters;
	}

	public void setFixedparameters(List fixedparameters) 
	{
		this.fixedParameters = fixedparameters;
	}

	public void addFixedParameter(FixedParameter fixedparameter) 
	{
		fixedParameters.add(fixedparameter);
		Collections.sort(fixedParameters);
	}

	public int compareTo(Object arg0) 
	{
		return _functionName.compareTo(((Function) arg0).get_functionName());
	}

	public String get_fortranFunction() {
		return _fortranFunction;
	}

	public void set_fortranFunction(String function) {
		_fortranFunction = function;
	}

	public boolean is_release() {
		return _release;
	}

	public void set_release(boolean _release) {
		this._release = _release;
	}

}

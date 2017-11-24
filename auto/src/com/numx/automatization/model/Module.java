package com.numx.automatization.model;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

public class Module implements Comparable {

	private int _moduleId;
	private String _moduleName;
	private String _moduleDescription ="";
	private boolean _onlyC;
	private boolean _release;
	private List functions;
	
	public Module(int id, String name, String description, boolean onlyc, boolean release) {
		_moduleId = id;
		_moduleName = name;
		_moduleDescription = description;
		_onlyC = onlyc;
		_release = release;
		functions = new ArrayList();
	}

	public int get_moduleId() {
		return _moduleId;
	}

	public void set_moduleId(int id) {
		_moduleId = id;
	}

	public String get_moduleDescription() {
		return _moduleDescription;
	}

	public void set_moduleDescription(String description) {
		_moduleDescription = description;
	}

	public String get_moduleName() {
		return _moduleName;
	}

	public String get_moduleNameNoSpace() {
		return _moduleName.replaceAll(" ","");
	}

	public String get_moduleNameLowerNoSpace() {
		return _moduleName.toLowerCase().replaceAll(" ","");
	}

	public String get_moduleNameLowerScore() {
		return _moduleName.toLowerCase().replaceAll(" ","-");
	}

	public String get_moduleNameLowerUnderScore() {
		return _moduleName.toLowerCase().replaceAll(" ","_");
	}

	public void set_moduleName(String name) {
		_moduleName = name;
	}

	public boolean is_onlyC() {
		return _onlyC;
	}

	public void set_onlyC(boolean _onlyc) {
		_onlyC = _onlyc;
	}
	
	/**
	 * Get functions ordered by name
	 * @return
	 */
	public List getFunctions() {
		return functions;
	}

	public void setFunctions(List functions) {
		this.functions = functions;
	}

	public void addFunction(Function function) {
		functions.add(function);
		Collections.sort(functions);
	}

	public int compareTo(Object arg0) {
		return _moduleName.compareTo(((Module) arg0).get_moduleName());	
	}

	public boolean is_release() {
		return _release;
	}

	public void set_release(boolean _release) {
		this._release = _release;
	}
}

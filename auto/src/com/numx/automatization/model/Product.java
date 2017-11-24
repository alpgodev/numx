package com.numx.automatization.model;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

public class Product implements Comparable {
	
	private int _productId;
	private String _productName;
	private String _productShortName;
	private String _productDescription = "";
	private String _productVersion;
	private List modules;

	public Product(int id, String name, String shortName, String description, String version) {
		_productId = id;
		_productName = name;
		_productShortName = shortName;
		_productDescription = description;
		_productVersion = version;
		modules = new ArrayList();
	}

	public String get_productDescription() {
		return _productDescription;
	}

	public void set_productDescription(String description) {
		_productDescription = description;
	}

	public String get_productName() {
		return _productName;
	}

	public void set_productName(String name) {
		_productName = name;
	}

	public String get_productShortName() {
		return _productShortName;
	}

	public String get_productShortNameLower() {
		return _productShortName.toLowerCase();
	}

	public void set_productShortName(String shortName) {
		_productShortName = shortName;
	}

	public String get_productVersion() {
		return _productVersion;
	}

	public void set_productVersion(String version) {
		_productVersion = version;
	}

	public int get_productId() {
		return _productId;
	}

	public void set_productId(int id) {
		_productId = id;
	}

	/**
	 * Gets modules ordered by product name
	 * @return
	 */
	public List getModules() {
		return modules;
	}

	public void setModules(List modules) {
		this.modules = modules;
	}

	public void addModule(Module module) {
		modules.add(module);
		Collections.sort(modules);
	}
	
	public int compareTo(Object arg0) {
		return this.get_productName().compareTo(((Product) arg0).get_productName());
	}
}
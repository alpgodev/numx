package com.numx.automatization.xml;

import java.io.File;
import java.io.PrintStream;
import java.util.Iterator;
import java.util.List;

import com.numx.automatization.model.Function;
import com.numx.automatization.model.Module;
import com.numx.automatization.model.Product;
import com.numx.automatization.model.Parameter;

public class XMLConstants {

	public final static String PRODUCTS_XML_FILE = "products.xml";

	public final static String PRODUCTS_NODE = "products";
	public final static String PRODUCT_NODE = "product";
	public final static String PRODUCT_NAME = "name";
	public final static String PRODUCT_SHORT_NAME = "short-name";
	public final static String PRODUCT_DESC = "description";
	public final static String PRODUCT_VERSION = "version";
		
	public final static String MODULE_NODE = "module";
	public final static String MODULE_NAME = "name";
	public final static String MODULE_C_ONLY = "c-only";
	public final static String MODULE_RELEASE = "release";
	public final static String MODULE_DESC = "description";

	public final static String FUNCTION_NODE = "function";
	public final static String FUNCTION_NAME = "name";
	public final static String FUNCTION_C_ONLY = "c-only";
	public final static String FUNCTION_RELEASE = "release";
	public final static String FUNCTION_DESC = "description";
	public final static String FUNCTION_FORTRAN = "fortran-function";
    public final static String FUNCTION_C = "c-function";
	public final static String FUNCTION_SHORT_DESC = "short-description";
	public final static String FUNCTION_IS_TEST = "is-test";
	public final static String FUNCTION_GW_FUNCTION_NAME = "gw-function-name";
	public final static String FUNCTION_IS_FORTRAN = "is-fortran";


	public final static String PARAMETER_NODE = "parameter";
	public final static String PARAMETER_NAME = "name";
	public final static String PARAMETER_TYPE = "type";
	public final static String PARAMETER_INPUT = "input-parameter";
	public final static String PARAMETER_OUTPUT_FUNCTION = "output-function";
	public final static String PARAMETER_LAST_INPUT = "last-input";
	public final static String PARAMETER_ROW = "row";
	public final static String PARAMETER_COLUMN = "column";
	public final static String PARAMETER_TEST_ROW = "testrow";
	public final static String PARAMETER_TEST_COLUMN = "testcol";
	public final static String PARAMETER_DESC = "description";
	public final static String PARAMETER_ORDER = "order";
	public final static String PARAMETER_ZERO_LENGTH = "zero-length";
	public final static String PARAMETER_RANDOM = "random";
	
	public final static String CONSTRAINT_NODE= "constraint";
	public final static String CONSTRAINT_TYPE = "type";
	public final static String CONSTRAINT_VALUE = "value";
	public final static String CONSTRAINT_VALUE2 = "value2";
	public final static String CONSTRAINT_VALUE_JAVA = "valueJava";
	public final static String CONSTRAINT_VALUE_JAVA2 = "valueJava2";
	public final static String CONSTRAINT_FIRST_INDEX = "firstIndex";
	public final static String CONSTRAINT_LAST_INDEX = "lastIndex";
	public final static String CONSTRAINT_FIRST_INDEX_JAVA = "firstIndexJava";
	public final static String CONSTRAINT_LAST_INDEX_JAVA = "lastIndexJava";

	public final static String SUM_NODE= "sum";
	public final static String SUM_TYPE = "type";
	public final static String SUM_VALUE = "value";
	public final static String SUM_NAME = "name";
	public final static String SUM_COND_VALUE1 = "condValue1";
	public final static String SUM_COND_VALUE2 = "condValue2";
	public final static String SUM_COND_TYPE = "condType";
	public final static String SUM_VALUE_JNI = "valueJNI";
	public final static String SUM_COND_VALUE1_JNI = "condValue1JNI";
	public final static String SUM_COND_VALUE2_JNI = "condValue2JNI";
	public final static String SUM_VALUE_VB = "valueVB";
	public final static String SUM_COND_VALUE1_VB = "condValue1VB";
	public final static String SUM_COND_VALUE2_VB = "condValue2VB";

	public final static String GETWORKSPACE_NODE= "gw";
	public final static String GETWORKSPACE_NAME = "name";
	public final static String GETWORKSPACE_TYPE = "type";
	public final static String GETWORKSPACE_PARAMETER_TYPE = "paramType";

	public final static String FIXEDPARAMETER_NODE= "fixedParameter";
	public final static String FIXEDPARAMETER_NAME = "name";
	public final static String FIXEDPARAMETER_INITIAL_VALUE = "value";
	public final static String FIXEDPARAMETER_TYPE = "type";

	public final static String FP_INITIALVALUE_NODE= "fpInitialValue";
	public final static String FP_INITIALVALUE_INDEX= "index";// If index <0 then the vector is initialized to the value="value"
	public final static String FP_INITIALVALUE_VALUE= "value";
	

	
	public static String getModuleFileName(Product product, Module module, File rootDir) {
		StringBuffer sb = new StringBuffer(rootDir.getAbsolutePath());
		sb.append('/');
		sb.append(product.get_productShortName());
		sb.append('-');
		sb.append(module.get_moduleName());
		sb.append(".xml");
		for(int i=0; i<sb.length(); i++) {
			if (Character.isSpaceChar(sb.charAt(i))) {
				sb.setCharAt(i, '_');
			}
		}
		return sb.toString();
	}

	public static void printProducts(PrintStream out, List products) {
		for (Iterator it = products.iterator(); it.hasNext();) {
			Product product = (Product) it.next();
			out.println("-Product: " + product.get_productName());
			for (Iterator it2 = product.getModules().iterator(); it2.hasNext();) {
				Module module = (Module) it2.next();
				out.println("\t-Module: " + module.get_moduleName());
				for (Iterator it3 = module.getFunctions().iterator(); it3.hasNext();) {
					Function f = (Function) it3.next();
					out.println("\t\t-function: " + f.get_functionName() 
							+ " (" + f.getParameters().size() +")" );
				}
			}
		}
	}

}

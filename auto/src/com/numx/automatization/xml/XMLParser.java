package com.numx.automatization.xml;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;

import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;
import org.w3c.dom.Text;


import com.numx.automatization.model.Constraint;
import com.numx.automatization.model.Sum;
import com.numx.automatization.model.Function;
import com.numx.automatization.model.GetWorkspace;
import com.numx.automatization.model.FixedParameter;
import com.numx.automatization.model.FPInitialValue;
import com.numx.automatization.model.Module;
import com.numx.automatization.model.Parameter;
import com.numx.automatization.model.Product;
import com.numx.automatization.xml.XMLConstants;

/**
 * Parse xml files for loading meta information that describes the product
 * suite. Meta information is stored in several xml files :
 * - a 'products.xml' file containing the product list and their module
 * - a xml file per module
 *  
 * @see XMLConstants
 * @author
 */
public class XMLParser extends XMLConstants {

	
	public static void main(String[] args) throws Exception {
        printProducts(System.out, new XMLParser().getProducts(args[0]));
	}

	public List getProducts(String outDir) throws Exception {
		return readProducts(new File(outDir));
	}
	
	/**
	 * Read the product meta data from XML files 
	 * @param fileName is the name of the main XML file
	 * @return a list of Product instances
	 * @throws Exception
	 */
	private List readProducts(File outDir) throws Exception {
		ArrayList products = new ArrayList();
        DocumentBuilderFactory domFactory = DocumentBuilderFactory.newInstance();
        
        //Read main XML document
        domFactory.setValidating(false);
        DocumentBuilder parser = domFactory.newDocumentBuilder();
       
        File productsFile = new File(outDir, PRODUCTS_XML_FILE);
		Document document = parser.parse(productsFile);
		
		Node productsNode = document.getFirstChild();
		Map productsChildren = groupChildrenByName(productsNode);
        List productNodes = (List) productsChildren.remove(PRODUCT_NODE);
        if (productNodes != null) {
	        for (int i = 0; i < productNodes.size(); i++) {
	            Node productNode = (Node) productNodes.get(i);
	            String name = getNodeAttribute(productNode, PRODUCT_NAME, null);
	            String shortName = getNodeAttribute(productNode, PRODUCT_SHORT_NAME, null);
	            String desc = getNodeAttribute(productNode, PRODUCT_DESC, "");
	            String version = getNodeAttribute(productNode, PRODUCT_VERSION, null);
	            Product product = new Product(i, name, shortName, desc, version);
	            
	    		Map productChildren = groupChildrenByName(productNode);
	            List moduleNodes = (List) productChildren.remove(MODULE_NODE);
	            if (moduleNodes != null) {
	    	        for (int j = 0; j < moduleNodes.size(); j++) {
	    	            Node moduleNode = (Node) moduleNodes.get(j);
	    	            String moduleName = getNodeAttribute(moduleNode, MODULE_NAME, null);
	    	            Module module = new Module(i*1000 + j, moduleName, null, false, true);
	    	            
	    	            //Read module info from dedicated xml file
	    	            readModule(product, module, outDir);
	    	            product.addModule(module);
	    	        }
	            }
	            products.add(product);
	        }
        }
        Collections.sort(products);
		return products;
	}


	/**
	 * Read the XML file corresponding to the module
	 * @param product is the product meta object
	 * @param module is the module meta object to fill
	 */
	private void readModule(Product product, Module module, File rootDir) throws Exception {
		String fileName = getModuleFileName(product, module, rootDir);
        DocumentBuilderFactory domFactory = DocumentBuilderFactory.newInstance();
        domFactory.setValidating(false);
        DocumentBuilder parser = domFactory.newDocumentBuilder();
		Document document = parser.parse(new File(fileName));
		Node moduleNode = document.getFirstChild();

		String desc = getSubTextNode(moduleNode, MODULE_DESC, "");
		module.set_moduleDescription(desc);
        
        String onlyC = getNodeAttribute(moduleNode, MODULE_C_ONLY, null);
        if (onlyC != null) {
        	module.set_onlyC(Boolean.valueOf(onlyC).booleanValue());
        }

        String release = getNodeAttribute(moduleNode, MODULE_RELEASE, null);
        if (release != null) {
        	module.set_release(Boolean.valueOf(release).booleanValue());
        }

		Map moduleChildren = groupChildrenByName(moduleNode);
        List functionNodes = (List) moduleChildren.remove(FUNCTION_NODE);
        if (functionNodes != null) {
	        for (int i = 0; i < functionNodes.size(); i++) {
	            module.addFunction(parseFunction((Node) functionNodes.get(i)));
	        }
        }
	}

	/**
	 * Parse DOM Node representing a function
	 * @param n is the DOM node
	 * @return a new Function meta object
	 */
	private Function parseFunction(Node n) {
		String name = getNodeAttribute(n, FUNCTION_NAME, null);
		String shortDesc = getSubTextNode(n, FUNCTION_SHORT_DESC, "");
		String desc = getSubTextNode(n, FUNCTION_DESC, "");
		String fortran = getNodeAttribute(n, FUNCTION_FORTRAN, null);
		String cOnly = getNodeAttribute(n, FUNCTION_C_ONLY, "false");
		String release = getNodeAttribute(n, FUNCTION_RELEASE, "true");
		String istest = getNodeAttribute(n, FUNCTION_IS_TEST, "false");
		String gwfunctionname = getNodeAttribute(n, FUNCTION_GW_FUNCTION_NAME, "");
		String isfortran = getNodeAttribute(n, FUNCTION_IS_FORTRAN, "true");

		Function f = new Function(0, name, shortDesc, desc, Boolean.valueOf(cOnly).booleanValue(), Boolean.valueOf(release).booleanValue(), 
			Boolean.valueOf(istest).booleanValue(), Boolean.valueOf(isfortran).booleanValue(), fortran, gwfunctionname);
		Map functionChildren = groupChildrenByName(n);
        List parameterNodes = (List) functionChildren.remove(PARAMETER_NODE);
        if (parameterNodes != null) {
	        for (int i = 0; i < parameterNodes.size(); i++) {
	            f.addParameter(parseParameter((Node) parameterNodes.get(i)));
	        }
        }
		List gwNodes = (List) functionChildren.remove(GETWORKSPACE_NODE);
		if (gwNodes != null) 
		{
			for (int i = 0; i < gwNodes.size(); i++) 
			{
				f.addGetWorkspace(parseGetWorkspace((Node) gwNodes.get(i)));
			}
		}
		
		List fixedParameterNodes = (List) functionChildren.remove(FIXEDPARAMETER_NODE);
		if (fixedParameterNodes != null) 
		{
			for (int i = 0; i < fixedParameterNodes.size(); i++) 
			{
				f.addFixedParameter(parseFixedParameter((Node) fixedParameterNodes.get(i)));
			}
		}

		return f;
	}
	
	/**
	 * Parse DOM Node representing a parameter
	 * @param n is the DOM node
	 * @return a new Parameter meta object
	 */
	private Parameter parseParameter(Node n) {
		String name = getNodeAttribute(n, PARAMETER_NAME, null);
		String type = getNodeAttribute(n, PARAMETER_TYPE, null);
		String in = getNodeAttribute(n, PARAMETER_INPUT, null);
		String lastInput = getNodeAttribute(n, PARAMETER_LAST_INPUT, "false");
		String outfunction = getNodeAttribute(n, PARAMETER_OUTPUT_FUNCTION, "false");
		String order = getNodeAttribute(n, PARAMETER_ORDER, null);
		String row = getNodeAttribute(n, PARAMETER_ROW, null);
		String col = getNodeAttribute(n, PARAMETER_COLUMN, null);
		String testrow = getNodeAttribute(n, PARAMETER_TEST_ROW, null);
		String testcol = getNodeAttribute(n, PARAMETER_TEST_COLUMN, null);
		String desc = getNodeAttribute(n, PARAMETER_DESC, "");
		String len = getNodeAttribute(n, PARAMETER_ZERO_LENGTH, "false");
		String ran = getNodeAttribute(n, PARAMETER_RANDOM, "false");
		Parameter p = new Parameter(0, name, type, desc, Boolean.valueOf(in).booleanValue(), 
				Boolean.valueOf(lastInput).booleanValue(), Boolean.valueOf(outfunction).booleanValue(), Integer.parseInt(order), row, col, testrow, testcol, 
				Boolean.valueOf(len).booleanValue(), Boolean.valueOf(ran).booleanValue());
		Map parameterChildren = groupChildrenByName(n);
		List constraintNodes = (List) parameterChildren.remove(CONSTRAINT_NODE);
		if (constraintNodes != null){
			for (int i = 0; i < constraintNodes.size(); i++){
				p.addConstraint(parseConstraint((Node) constraintNodes.get(i)));
			}
		}
		List sumNodes = (List) parameterChildren.remove(SUM_NODE);
		if (sumNodes != null)
		{
			for (int i = 0; i < sumNodes.size(); i++)
			{
				p.addSum(parseSum((Node) sumNodes.get(i)));
			}
		}
		return p;
	}
	
	/**
	 * Parse DOM Node representing a parameter
	 * @param n is the DOM node
	 * @return a new Parameter meta object
	 */
	private Sum parseSum(Node n) 
	{	
		String name = getNodeAttribute(n, SUM_NAME, null);
		String value = getNodeAttribute(n, SUM_VALUE, null);
		String valueJNI = getNodeAttribute(n, SUM_VALUE_JNI, null);
		String valueVB = getNodeAttribute(n, SUM_VALUE_VB, null);
		String type = getNodeAttribute(n, SUM_TYPE, null);
		String condValue1 = getNodeAttribute(n, SUM_COND_VALUE1, null);
		String condValue2 = getNodeAttribute(n, SUM_COND_VALUE2, null);
		String condValue1JNI = getNodeAttribute(n, SUM_COND_VALUE1_JNI, null);
		String condValue2JNI = getNodeAttribute(n, SUM_COND_VALUE2_JNI, null);
		String condValue1VB = getNodeAttribute(n, SUM_COND_VALUE1_VB, null);
		String condValue2VB = getNodeAttribute(n, SUM_COND_VALUE2_VB, null);
		String condType = getNodeAttribute(n, SUM_COND_TYPE, null);
		Sum s = new Sum(name, type, value, condValue1, condValue2, condType, 
						valueJNI, condValue1JNI, condValue2JNI,
						valueVB, condValue1VB, condValue2VB);
		return s;
	}


	/**
	 * Parse DOM Node representing a parameter
	 * @param n is the DOM node
	 * @return a new Parameter meta object
	 */
	private Constraint parseConstraint(Node n) 
	{
		String type = getNodeAttribute(n, CONSTRAINT_TYPE, null);
		String value = getNodeAttribute(n, CONSTRAINT_VALUE, null);
		String value2 = getNodeAttribute(n, CONSTRAINT_VALUE2, null);
		String valueJava = getNodeAttribute(n, CONSTRAINT_VALUE_JAVA, null);
		String valueJava2 = getNodeAttribute(n, CONSTRAINT_VALUE_JAVA2, null);
		String firstIndex = getNodeAttribute(n, CONSTRAINT_FIRST_INDEX, null);
		String lastIndex = getNodeAttribute(n, CONSTRAINT_LAST_INDEX, null);
		String firstIndexJava = getNodeAttribute(n, CONSTRAINT_FIRST_INDEX_JAVA, null);
		String lastIndexJava = getNodeAttribute(n, CONSTRAINT_LAST_INDEX_JAVA, null);
		Constraint c = new Constraint(type, value, value2, valueJava, valueJava2,
									firstIndex, lastIndex, firstIndexJava, lastIndexJava);
		return c;
	}

	/**
	 * Parse DOM Node representing a parameter
	 * @param n is the DOM node
	 * @return a new Parameter meta object
	 */
	private GetWorkspace parseGetWorkspace(Node n) 
	{
		String name = getNodeAttribute(n, GETWORKSPACE_NAME, null);
		String type = getNodeAttribute(n, GETWORKSPACE_TYPE, null);
		String paramType = getNodeAttribute(n, GETWORKSPACE_PARAMETER_TYPE, null);
		GetWorkspace gw = new GetWorkspace(name, type, paramType);
		return gw;
	}

	/**
	 * Parse DOM Node representing a parameter
	 * @param n is the DOM node
	 * @return a new Parameter meta object
	 */
	private FixedParameter parseFixedParameter(Node n) 
	{
		String name = getNodeAttribute(n, FIXEDPARAMETER_NAME, null);//
		String value = getNodeAttribute(n, FIXEDPARAMETER_INITIAL_VALUE, null);
		String type = getNodeAttribute(n, FIXEDPARAMETER_TYPE, null);
		FixedParameter fixedParameter = new FixedParameter(name, value, type);
		Map fixedParameterChildren = groupChildrenByName(n);
		List fpInitialValuesNodes = (List) fixedParameterChildren.remove(FP_INITIALVALUE_NODE);
		if (fpInitialValuesNodes != null)
		{
			for (int i = 0; i < fpInitialValuesNodes.size(); i++)
			{
				fixedParameter.addFPInitialValue(parseFPInitialValue((Node) fpInitialValuesNodes.get(i)));
			}
		}
		return fixedParameter;
	}

	/**
	 * Parse DOM Node representing a parameter
	 * @param n is the DOM node
	 * @return a new Parameter meta object
	 */
	private FPInitialValue parseFPInitialValue(Node n) 
	{
		
		String index = getNodeAttribute(n, FP_INITIALVALUE_INDEX, null);
		String value = getNodeAttribute(n, FP_INITIALVALUE_VALUE, null);
		FPInitialValue fpInitialValue = new FPInitialValue(value, index);
		return fpInitialValue;
	}

	/**
	 * Get the attribute value from a node
	 * @param n is the node hvaing attributes
	 * @param name is the attribute name
	 * @param defaultValue is the default value if the attribut does not exist.
	 * @return the string value found or the default value if not found
	 */
	private String getNodeAttribute(Node n, String name, String defaultValue) {
        Node n2 = n.getAttributes().getNamedItem(name);
        return (n2 == null ?  defaultValue:  n2.getNodeValue());
	}
	
	/**
	 * Get the Text value from a sub node
	 * @param n is the node hvaing a sub node
	 * @param name is the sub node name
	 * @param defaultValue is the default value if the sub node does not exist.
	 * @return the string value found or the default value if not found
	 */
	private String getSubTextNode(Node n, String name, String defaultValue) {
		Map children = groupChildrenByName(n);
		List l = (List) children.get(name);
		if (l == null || l.size() == 0) {
			return null;
		}
		if (l.size() > 1) {
			System.err.println("Several children '" + name 
					+ "' for the parent node '" + n .getNodeName() + "'.");
		}
		Text res = (Text)((Node) l.get(0)).getFirstChild();
		return res == null ? defaultValue : res.getData();
	}
	
	/**
	 * Sorts sub nodes by their name
	 * @param n is the node having sub nodes
	 * @return a map<String nodeName, List<Node sub nodes>>
	 */
    private Map groupChildrenByName(Node n) {
        NodeList nl = n.getChildNodes();
        int size = nl.getLength(); 
        if (size == 0) {
            return Collections.EMPTY_MAP;
        }
        HashMap result = new HashMap(size);
        for (int i = 0; i < size; i++) {
            Node child = nl.item(i);
            String name = child.getNodeName();
            List l = (List) result.get(name);
            if (l == null) {
                l = new ArrayList();
                result.put(name, l);
            }
            l.add(child);
        }
        return result;
    }

}

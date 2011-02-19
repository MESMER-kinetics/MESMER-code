
import java.awt.Component;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.event.ActionEvent;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.util.ArrayList;
import javax.swing.JButton;
import javax.swing.JComboBox;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.ScrollPaneConstants;

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

/**
 *
 * @author mhh
 */
public class MESMER_v3 extends JFrame {

    public String xMLStore = "";
    public int height = 200;
    public JPanel mainPanel = new JPanel();
    public JLabel mesmerLable = new JLabel();

    public JPanel panelStoreReaction = new JPanel();

    public JScrollPane scrollPanReactions = new JScrollPane();
    public JScrollPane scrollPanAll = new JScrollPane();

    public JLabel measuringUnitLabel = new JLabel();
    public JComboBox measuringUnitCombo = new JComboBox(new String[] {"kJ/mol", "kcal/mol", "hartree", "cm-1"} );

    public JButton addReaction = new JButton();
    public JButton viewAllRecords = new JButton();

    public int reactionNumber = 1;

    public JButton submitButton = new JButton();

    public static String measuringUnitValue = "kJ/mol";

    public ReactionContiontion reactionConditionObject = new ReactionContiontion();
    public  Controls controlObject = new Controls();

    public ViewAllRecordFrame viewAllRecordFrameObject;

    MESMER_v3()
    {
        
        
        setSize(850,700);
        setDefaultCloseOperation ( JFrame.EXIT_ON_CLOSE );
        setLocationRelativeTo(null);
        setTitle("MESMER");
        setResizable(false);
        setVisible(true);


        //Main panle

       mainPanel.setLayout(new FlowLayout(FlowLayout.CENTER ));
       mainPanel.setSize(830, 1000);
       mainPanel.setPreferredSize(new Dimension(830, 1000));

       // mesmer lable
       mesmerLable.setFont(new java.awt.Font("Wide Latin", 0, 18));
      // mesmerLable.setHorizontalAlignment(SwingConstants.RIGHT);
       // this is to space the title in the left.  May be better to set the flow layout
       // manager to the left?
       String space = "                                                     ";
       mesmerLable.setText(space + space + "MESMER");
//         mesmerLable.setText("MESMER");
       mesmerLable.setSize(807, 40);
       mesmerLable.setPreferredSize(new Dimension(807,60));
      // mesmerLable.setPreferredSize(new Dimension(807, 60));

       mainPanel.add(mesmerLable);
       mainPanel.setVisible(true);

       //set the panle to store the reaction data
       
       panelStoreReaction.setLayout(new FlowLayout());
       panelStoreReaction.setPreferredSize(new Dimension(730,180));
       panelStoreReaction.add(new ReactionPan2( 1));
      // panelStoreReaction.add(new ReactionPan2( 2));

      scrollPanReactions.setHorizontalScrollBarPolicy(ScrollPaneConstants.HORIZONTAL_SCROLLBAR_ALWAYS);
      scrollPanReactions.setVerticalScrollBarPolicy(ScrollPaneConstants.VERTICAL_SCROLLBAR_ALWAYS);
      scrollPanReactions.setViewportView(panelStoreReaction);
         //  scrollPan.add(mainJPane2);
      scrollPanReactions.setPreferredSize(new Dimension(780,200));
      mainPanel.add(scrollPanReactions);
    

       measuringUnitLabel.setFont(new java.awt.Font("Tahoma", 1, 12)); // NOI18N
       measuringUnitLabel.setText("Measuring Units");
       mainPanel.add(measuringUnitLabel);
 //measuringUnitLabel.setBounds(60, 260, 110, 30);
       measuringUnitCombo.setSize(61, 20);
       measuringUnitCombo.setPreferredSize(new Dimension(111,20));
       
        measuringUnitCombo.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed (java.awt.event.ActionEvent evt) {
                measuringUnitComboActionPerformed(evt);
            }
        });
       
       mainPanel.add(measuringUnitCombo);
       
       

       addReaction.setText("Adfd New Reaction");
       viewAllRecords.setText("View All Records");

       addReaction.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                addReactionActionPerformed(evt);

            }
        });

        viewAllRecords.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                viewAllRecordsActionPerformed(evt);

            }
        });     
        
       mainPanel.add(addReaction);
       mainPanel.add(viewAllRecords);
       mainPanel.add(reactionConditionObject);
       //mainPanel.add(new Controls());
      mainPanel.add(controlObject);


       submitButton.setText("Submit Job");

       submitButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                submitButtonsActionPerformed(evt);

            }
        });

       submitButton.setSize(300, 20);
       submitButton.setPreferredSize(new Dimension(300, 20));
       mainPanel.add(submitButton);

      scrollPanAll.setHorizontalScrollBarPolicy(ScrollPaneConstants.HORIZONTAL_SCROLLBAR_ALWAYS);
      scrollPanAll.setVerticalScrollBarPolicy(ScrollPaneConstants.VERTICAL_SCROLLBAR_ALWAYS);
      scrollPanAll.setViewportView(mainPanel);
         //  scrollPan.add(mainJPane2);
      scrollPanAll.setPreferredSize(new Dimension(830,500));
    //  mainPanel.add(scrollPanReactions);

       
        add(scrollPanAll);
    }


    public void addReactionActionPerformed(ActionEvent evt)
    {
        boolean add = false;
        int counter, componentsToCount;
        String errorOccuredAt ="";
        counter = 0;
        Component[] components = panelStoreReaction.getComponents();
        componentsToCount = components.length;
        System.out.println("The componentCounter is " + componentsToCount );
        ReactionPan2 reactantValues;
        for(int i =0; i< components.length; i++)
        {
            reactantValues = (ReactionPan2) components[i];
            if(reactantValues.getNo_ReactantComboSelected().equals("2") &&
                    reactantValues.getReactantCombo1().equals("Deficent") &&
                     reactantValues.getReactantCombo2().equals("Excess") &&
                     reactantValues.getNo_ProductComboSelected().equals("2") &&
                     reactantValues.getProductCombo1().equals("Deficent") &&
                     reactantValues.getProductCombo2().equals("Excess") &&

                     reactantValues.get_RRKM() == true)
            {
                counter++;
                System.out.println("Reactant 2, deficient/ ecess.  Product 2, deficient/excess."+
                        "And RRKM");
            }

            else if(reactantValues.getNo_ReactantComboSelected().equals("2") &&
            reactantValues.getReactantCombo1().equals("Deficent") &&
             reactantValues.getReactantCombo2().equals("Excess") &&
             reactantValues.getNo_ProductComboSelected().equals("1") &&
             reactantValues.getProductCombo1().equals("Model") &&

             reactantValues.get_RRKM() == false)
            {
                counter++;
                System.out.println("Reactant 2, deficient/ ecess.  Product 1, excess."+
                        "And ILT");
            }

            else if(reactantValues.getNo_ReactantComboSelected().equals("1") &&
                    reactantValues.getReactantCombo1().equals("Model") &&
                     reactantValues.getNo_ProductComboSelected().equals("2") &&
                     reactantValues.getProductCombo1().equals("Sink") &&
                     reactantValues.getProductCombo2().equals("Sink") &&

                     reactantValues.get_RRKM() == false)
            {
                counter++;
                System.out.println("Reactant 2, deficient/ ecess.  Product 2, deficient/excess."+
                        "And RRKM");
            }

            else if(reactantValues.getNo_ReactantComboSelected().equals("1") &&
                    reactantValues.getReactantCombo1().equals("Model") &&                     
                     reactantValues.getNo_ProductComboSelected().equals("1") &&
                     reactantValues.getProductCombo1().equals("Model") &&


                     reactantValues.get_RRKM() == true)
            {
                counter++;
                System.out.println("Reactant 2, deficient/ ecess.  Product 2, deficient/excess."+
                        "And RRKM");
            }

            else if(reactantValues.getNo_ReactantComboSelected().equals("2") &&
                    reactantValues.getReactantCombo1().equals("Model") &&
                     reactantValues.getReactantCombo2().equals("Model") &&
                     reactantValues.getNo_ProductComboSelected().equals("1") &&
                     reactantValues.getProductCombo1().equals("Model") &&


                     reactantValues.get_RRKM() == true)
            {
                counter++;
                System.out.println("Reactant 2, deficient/ ecess.  Product 2, deficient/excess."+
                        "And RRKM");
            }

            else if(reactantValues.getNo_ReactantComboSelected().equals("1") &&
            reactantValues.getReactantCombo1().equals("Model") &&
            reactantValues.getNo_ProductComboSelected().equals("2") &&
            reactantValues.getProductCombo1().equals("Model") &&
            reactantValues.getProductCombo1().equals("Model") &&

             reactantValues.get_RRKM() == true)
            {
                counter++;
                System.out.println("Reactant 2, deficient/ ecess.  Product 2, deficient/excess."+
                        "And RRKM");
            }

            else if(reactantValues.getNo_ReactantComboSelected().equals("1") &&
            reactantValues.getReactantCombo1().equals("Model") &&
            reactantValues.getNo_ProductComboSelected().equals("2") &&
            reactantValues.getProductCombo1().equals("Sink") &&
            reactantValues.getProductCombo1().equals("Sink") &&

             reactantValues.get_RRKM() == true)
            {
                counter++;
                System.out.println("Reactant 2, deficient/ ecess.  Product 2, deficient/excess."+
                        "And RRKM");
            }

            else {
                String addAtTheEnds = (i == (components.length -1))? " ": ", ";
                errorOccuredAt += String.valueOf((i+1)) + addAtTheEnds;

            }


        }

        if (counter == componentsToCount)
        {
            add = true;
        }
        else
        {
            JOptionPane.showMessageDialog(null,
            "You have an error in your reaction combination at locations " + errorOccuredAt ,
            "Reaction Combination Incorrect",
            JOptionPane.WARNING_MESSAGE);
        }


        if(add == true)
        {
             height += 190;
             reactionNumber++;
             panelStoreReaction.add(new  ReactionPan2(reactionNumber));
             panelStoreReaction.setPreferredSize(new Dimension(740,height));
             scrollPanReactions.setViewportView(panelStoreReaction);
        }
    }

    public void measuringUnitComboActionPerformed(ActionEvent evt)
    {
      measuringUnitValue = (String)measuringUnitCombo.getSelectedItem();
      System.out.println("THe measuring unit selected is " + measuringUnitValue);
    }

    public void  viewAllRecordsActionPerformed(ActionEvent evt)
    {
            Component[] components = panelStoreReaction.getComponents();
            Controls gasBufferValueObejct = controlObject;
            String gasBufferSelectedValue = gasBufferValueObejct.get_BufferGas();
            System.out.println("Before parsing " + gasBufferSelectedValue );
        //    gasBufferSelectedValue =  gasBufferSelectedValue.equals("He - Helium")? "He" : "N2";
            System.out.println("After parsing " + gasBufferSelectedValue );
            viewAllRecordFrameObject = new ViewAllRecordFrame(components, measuringUnitValue, gasBufferSelectedValue );
            viewAllRecordFrameObject.setVisible(true);
    }

    public void  submitButtonsActionPerformed(ActionEvent evt)
    {

    Component[] components = panelStoreReaction.getComponents();

    ReactionPan2 reactantValues;
    for(int i =0; i< components.length; i++)
    {
       reactantValues = (ReactionPan2) components[i];

       if(reactantValues.getNo_ReactantComboSelected().equals("1") )
       {
           System.out.println("############################# Reactionm1 ###############################################");
           System.out.println("Reaction ");
           System.out.println( "************** " + i + "******************");
           System.out.println("******Reactant value is ********" + reactantValues.getReactant_1TextField());
           System.out.println( "******Energy ***** " +  i + " "+  reactantValues.reactant1Input.getGetEnergy());
           System.out.println("******Rotational Constants " + reactantValues.reactant1Input.getRotationalConstants());
           System.out.println( "****** Vibrational Frequency " +  reactantValues.reactant1Input.getVibrationalFrequecy());
           System.out.println( "******Scalar factor " +  reactantValues.reactant1Input.getScalingFactor() );
           System.out.println( "****** Symetry" +  reactantValues.reactant1Input.getSymmerty() );
           System.out.println( "******Moledule Weight " +  reactantValues.reactant1Input.getMolecularWeight() );
           System.out.println( "******Spin Multiplicity " +  reactantValues.reactant1Input.getSpinMulticity() );
           System.out.println(  "******Density of State " +  reactantValues.reactant1Input.getDensityOfState() );
           System.out.println( "******Ession " +  reactantValues.reactant1Input.getEpsion() );
           System.out.println( "******Sigma " +  reactantValues.reactant1Input.getSigma() );
           System.out.println( "****** Down" +  reactantValues.reactant1Input.getDown() );
           System.out.println( "****** Scientifi Rationa" +  reactantValues.reactant1Input.getScienfifiRational() );
       }

        if(reactantValues.getNo_ReactantComboSelected().equals("2") )
       {
            System.out.println("############################## Reaction 2 ##############################################");
            System.out.println("Reaction2 ");

            System.out.println("****** Reactant 1 " + reactantValues.getReactant_1TextField());
            System.out.println("****** Combo Selected 1 " + reactantValues.getReactantCombo1());

            System.out.println("****** Energy" + reactantValues.reactant2Input.getEnergy_1() );
            System.out.println("****** RotaitonalConstant " + reactantValues.reactant2Input.getRotationalConstant_1() );
            System.out.println("****** Vib Frequency " + reactantValues.reactant2Input.getVibrationalFrequency_1() );
            System.out.println("****** Scaling Factor" + reactantValues.reactant2Input.getScalingFactor_1() );
            System.out.println("****** Symertry " + reactantValues.reactant2Input.getSymmetry_1() );
            System.out.println("****** Moleculeare weight " + reactantValues.reactant2Input.getMolecularWeight_1() );
            System.out.println("******  Multiplicity " + reactantValues.reactant2Input.getSpinMultiplicity_1() );
            System.out.println("****** Density of State " + reactantValues.reactant2Input.getDensityofState_1() );
            System.out.println("******Scientific Ractional " + reactantValues.reactant2Input.getScientificRational_1() );

            System.out.println("****** Reactant 2 " + reactantValues.getReactant_2TextField());
            System.out.println("****** Combo Selected 1 " + reactantValues.getReactantCombo2());
            System.out.println("****** Energy " +  reactantValues.reactant2Input.getEnergy_2() );
            System.out.println("****** RotaitonalConstant " + reactantValues.reactant2Input.getRotationalConstant_2() );
            System.out.println("******Vib Frequency " + reactantValues.reactant2Input.getVibrationalFrequency_2() );
            System.out.println("****** Scaling Factor " + reactantValues.reactant2Input.getScalingFactor_2() );
            System.out.println("****** Symertry " + reactantValues.reactant2Input.getSymmetry_2() );
            System.out.println("****** Moleculeare weight  " + reactantValues.reactant2Input.getMolecularWeight_2() );
            System.out.println("****** Multiplicity " + reactantValues.reactant2Input.getSpinMultiplicity_2() );
            System.out.println("****** Density of State " + reactantValues.reactant2Input.getDensityofState_2() );
            System.out.println("****** Scientific Ractiona " + reactantValues.reactant2Input.getScientificRational_2() );
        }

       if(reactantValues.getNo_ProductComboSelected().equals("1"))
       {
            System.out.println("############################Product 1 ################################################");
            System.out.println("Product ");

            System.out.println( "************** " + i + "******************");
            System.out.println("******Reactant value is ********" + reactantValues.getProduct_1TextField());
            System.out.println( "******Energy ***** " +  i + " "+  reactantValues.product1Input.getGetEnergy());
            System.out.println("******Rotational Constants " + reactantValues.product1Input.getRotationalConstants());
            System.out.println( "****** Vibrational Frequency " +  reactantValues.product1Input.getVibrationalFrequecy());
            System.out.println( "******Scalar factor " +  reactantValues.product1Input.getScalingFactor() );
            System.out.println( "****** Symetry" +  reactantValues.product1Input.getSymmerty() );
            System.out.println( "******Moledule Weight " +  reactantValues.product1Input.getMolecularWeight() );
            System.out.println( "******Spin Multiplicity " +  reactantValues.product1Input.getSpinMulticity() );
            System.out.println(  "******Density of State " +  reactantValues.product1Input.getDensityOfState() );
            System.out.println( "******Ession " +  reactantValues.product1Input.getEpsion() );
            System.out.println( "******Sigma " +  reactantValues.product1Input.getSigma() );
            System.out.println( "****** Down" +  reactantValues.product1Input.getDown() );
            System.out.println( "****** Scientifi Rationa" +  reactantValues.product1Input.getScienfifiRational() );
       }

        if(reactantValues.getNo_ProductComboSelected().equals("2"))
       {
            System.out.println("############################Product 2 ################################################");
            System.out.println("Product 2 ");

            System.out.println("****** Product 1 " + reactantValues.getProduct_1TextField());
            System.out.println(" ***** Combo 1 selected " +  reactantValues.getProductCombo1());

            System.out.println("****** Energy" + reactantValues.product2Input.getEnergy_1() );
            System.out.println("****** RotaitonalConstant " + reactantValues.product2Input.getRotationalConstant_1() );
            System.out.println("****** Vib Frequency " + reactantValues.product2Input.getVibrationalFrequency_1() );
            System.out.println("****** Scaling Factor" + reactantValues.product2Input.getScalingFactor_1() );
            System.out.println("****** Symertry " + reactantValues.product2Input.getSymmetry_1() );
            System.out.println("****** Moleculeare weight " + reactantValues.product2Input.getMolecularWeight_1() );
            System.out.println("******  Multiplicity " + reactantValues.product2Input.getSpinMultiplicity_1() );
            System.out.println("****** Density of State " + reactantValues.product2Input.getDensityofState_1() );
            System.out.println("******Scientific Ractional " + reactantValues.product2Input.getScientificRational_1() );

            System.out.println("****** Product 2 " + reactantValues.getProduct_2TextField());
            System.out.println(" ***** Combo 2 selected " +  reactantValues.getProductCombo2());
            System.out.println("****** Energy " +  reactantValues.product2Input.getEnergy_2() );
            System.out.println("****** RotaitonalConstant " + reactantValues.product2Input.getRotationalConstant_2() );
            System.out.println("******Vib Frequency " + reactantValues.product2Input.getVibrationalFrequency_2() );
            System.out.println("****** Scaling Factor " + reactantValues.product2Input.getScalingFactor_2() );
            System.out.println("****** Symertry " + reactantValues.product2Input.getSymmetry_2() );
            System.out.println("****** Moleculeare weight  " + reactantValues.product2Input.getMolecularWeight_2() );
            System.out.println("****** Multiplicity " + reactantValues.product2Input.getSpinMultiplicity_2() );
            System.out.println("****** Density of State " + reactantValues.product2Input.getDensityofState_2() );
            System.out.println("****** Scientific Ractiona " + reactantValues.product2Input.getScientificRational_2() );
        }

       if(reactantValues.get_RRKM() == false)
       {
           System.out.println("########################### ILT #################################################");
           System.out.println("ILT ");

           System.out.println("****** ILT A ********" + reactantValues.iltMethodInput.get_ILT_A());
           System.out.println("******  ILT N********" + reactantValues.iltMethodInput.get_ILT_N());
           System.out.println("****** ILT E ********" + reactantValues.iltMethodInput.get_ILT_E());
           System.out.println("****** Execc Resional COncentrationm ********" + reactantValues.iltMethodInput.get_ExcessRegionalConcentration());
        //   System.out.println("****** Excess ********" + reactantValues.iltMethodInput.get_Excess());
           System.out.println("****** Scientific Rationale ********" + reactantValues.iltMethodInput.get_ScienfifiRational());
       }

      if(reactantValues.get_RRKM() == true)
       {
           System.out.println("########################### Transition #################################################");
           System.out.println("Transition  ");
           System.out.println("****** " + reactantValues.get_TransitionStateValue() );
           System.out.println("****** " + reactantValues.transitionStateInput.getGetEnergy() );
           System.out.println("****** " + reactantValues.transitionStateInput.getRotationalConstants());
           System.out.println("****** " + reactantValues.transitionStateInput.getVibrationalFrequecy() );
           System.out.println("****** " + reactantValues.transitionStateInput.getScalingFactor() );
           System.out.println("****** " + reactantValues.transitionStateInput.getSymmerty() );
           System.out.println("****** " + reactantValues.transitionStateInput.getMolecularWeight() );
           System.out.println("****** " + reactantValues.transitionStateInput.getSpinMulticity() );
           System.out.println("****** " + reactantValues.transitionStateInput.getDensityOfState() );
           System.out.println("****** " + reactantValues.transitionStateInput.getImaginitiveFrequency() );
           System.out.println("****** " + reactantValues.transitionStateInput.getScienfifiRational() );
      }
    }


    //Nopw display the ReactionConditions
    System.out.println("************** ReactionConditions********************");
    ArrayList<String> pressure = reactionConditionObject.get_PressureValues();
    ArrayList<String> temperature = reactionConditionObject.get_TemperatureValues();
    System.out.println("Pressure values are " + pressure );
    System.out.println("Temperature values are " +  temperature);
    System.out.println("Units " + reactionConditionObject.get_UnitCalculation());
    System.out.println("Precision " + reactionConditionObject.get_Precision() );


     System.out.println("************** Controlls ********************");
     System.out.println(" get_calculateRateCoefficientsOnly " + controlObject.get_calculateRateCoefficientsOnly()  );
     System.out.println(" get_printReactionOperatorColumnSums " + controlObject.get_printReactionOperatorColumnSums()  );
     System.out.println(" get_printGrainkbE " + controlObject.get_printGrainkbE()  );
     System.out.println(" get_printReactionOperatorSize " + controlObject.get_printReactionOperatorSize()  );
     System.out.println(" get_printTunnellingCoefficients " + controlObject.get_printTunnellingCoefficients()  );
     System.out.println(" get_testRateConstants " + controlObject.get_testRateConstants()  );
     System.out.println(" get_printCellDOS " + controlObject.get_printCellDOS()  );
     System.out.println(" get_printGrainBoltzmann " + controlObject.get_printGrainBoltzmann()  );
     System.out.println(" get_printGrainkfE " + controlObject.get_printGrainkfE()  );
     System.out.println(" get_printSpeciesProfile  " + controlObject.get_printSpeciesProfile()   );
     System.out.println(" get_testDOS " + controlObject.get_testDOS()  );
     System.out.println(" get_useTheSameCellNumberForAllConditions " + controlObject.get_useTheSameCellNumberForAllConditions()  );
     System.out.println(" get_printCellTransitionStateFlux " + controlObject.get_printCellTransitionStateFlux()  );
     System.out.println(" get_printGrainDOS " + controlObject.get_printGrainDOS()  );
     System.out.println(" get_printGrainTransitionStateFlux " + controlObject.get_printGrainTransitionStateFlux()  );
     System.out.println(" get_printTunnelingCoefficients " + controlObject.get_printTunnelingCoefficients()  );
     System.out.println(" get_testMicroRates " + controlObject.get_testMicroRates()  );
     System.out.println("get_Eigenvalues  " + controlObject.get_Eigenvalues());
     System.out.println("get_BufferGas  " + controlObject.get_BufferGas());
     System.out.println("get_CalcMethod  " + controlObject.get_CalcMethod());

        //Now display the controls values


   //   Component allReaction = new panelStoreReaction.getComponents();



     //THis to check if the Eligevalue is an int and other values;
     //this will be false if an excelption is throw
       boolean proceed = true;

       String eigenvaluesString = controlObject.getEigenvalues();

       try {
            int eigenvaluesInteger = Integer.parseInt(eigenvaluesString);
             System.out.println(eigenvaluesInteger);
       }
       catch(NumberFormatException nFE) {
        System.out.println("Not an Integer for Eigenvalues " + nFE);
        JOptionPane.showMessageDialog(null,
            "The Eigen Values must be of type integer or the system will not run",
            "Eigen Values",
            JOptionPane.WARNING_MESSAGE);

            //an exceptiuon has been throw so change to false
            proceed = false;
    }

       if(proceed == true)
       {
            generateXML();
       }
 }

     public  String convertRecProcType(String value)
        {
            if(value.equals("Deficent"))
               return "deficentReactant";
            else if(value.equals("Excess"))
               return "excessReactant";
            else if(value.equals("Sink"))
              return "sink";
            else if(value.equals("Model"))
              return "modelled";
            else
              return "reactant";
        }


    public void generateXML()
    {

        String header = "<?xml version=\"1.0\" encoding=\"utf-8\" ?>" + "\n" +
        "<?xml-stylesheet type=\'text/xsl\' href=\'../../mesmer2.xsl\' media=\'other\'?> " + "\n" +
        "<?xml-stylesheet type=\'text/xsl\' href=\'../../mesmer1.xsl\' media='screen\'?> " + "\n" +
        "<me:mesmer xmlns=\"http://www.xml-cml.org/schema\" xmlns:me=\"http://www.chem.leeds.ac.uk/mesmer\">" + "\n" ;

//        StringBuffer xMLStore = new StringBuffer();
//
//        xMLStore.append(header);


        Component[] components = panelStoreReaction.getComponents();
        Component[] components2 = panelStoreReaction.getComponents();

        ReactionPan2 reactantValues;

        String moleculeValue; //stores all the molecules values for XML
        moleculeValue="";
        moleculeValue +=  "<title>" + "TITLE"  + "</title>" + "\n";
        moleculeValue +=  "  <moleculeList>"+ "\n";
        for(int i =0; i< components.length; i++)
        {
            reactantValues = (ReactionPan2) components[i];

            if(reactantValues.getNo_ReactantComboSelected().equals("1") )
            {

                moleculeValue += "   <molecule id=\""+ reactantValues.getReactant_1TextField() + "\">" + "\n" +
                "     <propertyList>"+ "\n" +
                "       <property dictRef=\"me:ZPE\">"+ "\n" +
                "          <scalar units=\""+ measuringUnitValue +"\">" + reactantValues.reactant1Input.getGetEnergy() + "</scalar>"+ "\n" +
                "           </property>"+ "\n" +
                        "<property dictRef=\"me:rotConsts\"> " + "\n" +
                          "<array units=\"cm-1\">" + reactantValues.reactant1Input.getRotationalConstants() +"</array> " + "\n" +
                        "</property> "+ "\n" +
                        "<property dictRef=\"me:vibFreqs\"> " + "\n" +
                          "<array units=\"cm-1\">" + reactantValues.reactant1Input.getVibrationalFrequecy() + "</array> " + "\n" +
                        "</property> "+ "\n" +
                       "<property dictRef=\"me:frequenciesScaleFactor\">" + "\n"+
                        "<scalar>"+   reactantValues.reactant1Input.getScalingFactor()+  "</scalar>" + "\n"+
                       "</property>" + "\n"+
                        "<property dictRef=\"me:symmetryNumber\">" + "\n" +
                          "<scalar>" + reactantValues.reactant1Input.getSymmerty()+ "</scalar> "+ "\n" +
                        "</property> " + "\n" +                        
                        "<property dictRef=\"me:MW\"> " + "\n" +
                          "<scalar units=\"amu\">" + reactantValues.reactant1Input.getMolecularWeight() + "</scalar> " + "\n" +
                        "</property> " + "\n" +
                         "<property dictRef=\"me:spinMultiplicity\"> " + "\n" +
                          "<scalar>" + reactantValues.reactant1Input.getSpinMulticity() + "</scalar> " + "\n" +
                        "</property> " + "\n" +
                        "<property dictRef=\"me:epsilon\"> " + "\n" +
                          "<scalar>" + reactantValues.reactant1Input.getEpsion() + "</scalar> " + "\n" +
                        "</property> " + "\n" +
                        "<property dictRef=\"me:sigma\"> " + "\n" +
                          "<scalar>" +reactantValues.reactant1Input.getSigma() + "</scalar>  " + "\n" +
                        "</property> " + "\n" +
                        "<property dictRef=\"me:deltaEDown\"> " + "\n" +
                          "<scalar units=\"cm-1\">" +  reactantValues.reactant1Input.getDown()+ " </scalar> " + "\n" +
                        "</property> " + "\n" +                       
                      "</propertyList> " + "\n" +
                      "<me:DOSCMethod>" + reactantValues.reactant1Input.getDensityOfState()+ "</me:DOSCMethod> " + "\n" +
                    "</molecule> " + "\n";
            }

              if(reactantValues.getNo_ReactantComboSelected().equals("2") )
               {
                     moleculeValue += "   <molecule id=\""+ reactantValues.getReactant_1TextField() + "\">" + "\n" +
                "     <propertyList>"+ "\n" +
                "       <property dictRef=\"me:ZPE\">"+ "\n" +
                "          <scalar units=\""+ measuringUnitValue +"\">" + reactantValues.reactant2Input.getEnergy_1() + "</scalar>"+ "\n" +
                "           </property>"+ "\n" +
                        "<property dictRef=\"me:rotConsts\"> " + "\n" +
                          "<array units=\"cm-1\">" + reactantValues.reactant2Input.getRotationalConstant_1()+"</array> " + "\n" +
                        "</property> "+ "\n" +
                        "<property dictRef=\"me:vibFreqs\"> " + "\n" +
                          "<array units=\"cm-1\">" + reactantValues.reactant2Input.getVibrationalFrequency_1() + "</array> " + "\n" +
                        "</property> "+ "\n" +
                       "<property dictRef=\"me:frequenciesScaleFactor\">" + "\n"+
                        "<scalar>"+   reactantValues.reactant2Input.getScalingFactor_1() +  "</scalar>" + "\n"+
                       "</property>" + "\n"+
                        "<property dictRef=\"me:symmetryNumber\">" + "\n" +
                          "<scalar>" + reactantValues.reactant2Input.getSymmetry_1() + "</scalar> "+ "\n" +
                        "</property> " + "\n" +
                        "<property dictRef=\"me:MW\"> " + "\n" +
                          "<scalar units=\"amu\">" + reactantValues.reactant2Input.getMolecularWeight_1() + "</scalar> " + "\n" +
                        "</property> " + "\n" +
                         "<property dictRef=\"me:spinMultiplicity\"> " + "\n" +
                          "<scalar>" + reactantValues.reactant2Input.getSpinMultiplicity_1() + "</scalar> " + "\n" +
                        "</property> " + "\n" +
                      "</propertyList> " + "\n" +
                      "<me:DOSCMethod>" + reactantValues.reactant2Input.getDensityofState_1()+ "</me:DOSCMethod> " + "\n" +
                    "</molecule> " + "\n";

                     // the second recatnat
                    moleculeValue += "   <molecule id=\""+ reactantValues.getReactant_2TextField() + "\">" + "\n" +
                "     <propertyList>"+ "\n" +
                "       <property dictRef=\"me:ZPE\">"+ "\n" +
                "          <scalar units=\""+ measuringUnitValue +"\">" + reactantValues.reactant2Input.getEnergy_2() + "</scalar>"+ "\n" +
                "           </property>"+ "\n" +
                        "<property dictRef=\"me:rotConsts\"> " + "\n" +
                          "<array units=\"cm-1\">" + reactantValues.reactant2Input.getRotationalConstant_2()+"</array> " + "\n" +
                        "</property> "+ "\n" +
                        "<property dictRef=\"me:vibFreqs\"> " + "\n" +
                          "<array units=\"cm-1\">" + reactantValues.reactant2Input.getVibrationalFrequency_2() + "</array> " + "\n" +
                        "</property> "+ "\n" +
                       "<property dictRef=\"me:frequenciesScaleFactor\">" + "\n"+
                        "<scalar>"+   reactantValues.reactant2Input.getScalingFactor_2() +  "</scalar>" + "\n"+
                       "</property>" + "\n"+
                        "<property dictRef=\"me:symmetryNumber\">" + "\n" +
                          "<scalar>" + reactantValues.reactant2Input.getSymmetry_2() + "</scalar> "+ "\n" +
                        "</property> " + "\n" +
                        "<property dictRef=\"me:MW\"> " + "\n" +
                          "<scalar units=\"amu\">" + reactantValues.reactant2Input.getMolecularWeight_2() + "</scalar> " + "\n" +
                        "</property> " + "\n" +
                         "<property dictRef=\"me:spinMultiplicity\"> " + "\n" +
                          "<scalar>" + reactantValues.reactant2Input.getSpinMultiplicity_2() + "</scalar> " + "\n" +
                        "</property> " + "\n" +
                      "</propertyList> " + "\n" +
                      "<me:DOSCMethod>" + reactantValues.reactant2Input.getDensityofState_2()+ "</me:DOSCMethod> " + "\n" +
                    "</molecule> " + "\n";
              }

           if(reactantValues.getNo_ProductComboSelected().equals("1"))
            {
                moleculeValue += "   <molecule id=\""+ reactantValues.getProduct_1TextField() + "\">" + "\n" +
                "     <propertyList>"+ "\n" +
                "       <property dictRef=\"me:ZPE\">"+ "\n" +
                "          <scalar units=\""+ measuringUnitValue +"\">" + reactantValues.product1Input.getGetEnergy() + "</scalar>"+ "\n" +
                "           </property>"+ "\n" +
                        "<property dictRef=\"me:rotConsts\"> " + "\n" +
                          "<array units=\"cm-1\">" + reactantValues.product1Input.getRotationalConstants() +"</array> " + "\n" +
                        "</property> "+ "\n" +
                        "<property dictRef=\"me:vibFreqs\"> " + "\n" +
                          "<array units=\"cm-1\">" + reactantValues.product1Input.getVibrationalFrequecy() + "</array> " + "\n" +
                        "</property> "+ "\n" +
                       "<property dictRef=\"me:frequenciesScaleFactor\">" + "\n"+
                        "<scalar>"+   reactantValues.product1Input.getScalingFactor()+  "</scalar>" + "\n"+
                       "</property>" + "\n"+
                        "<property dictRef=\"me:symmetryNumber\">" + "\n" +
                          "<scalar>" + reactantValues.product1Input.getSymmerty()+ "</scalar> "+ "\n" +
                        "</property> " + "\n" +
                        "<property dictRef=\"me:MW\"> " + "\n" +
                          "<scalar units=\"amu\">" + reactantValues.product1Input.getMolecularWeight() + "</scalar> " + "\n" +
                        "</property> " + "\n" +
                         "<property dictRef=\"me:spinMultiplicity\"> " + "\n" +
                          "<scalar>" + reactantValues.product1Input.getSpinMulticity() + "</scalar> " + "\n" +
                        "</property> " + "\n" +
                        "<property dictRef=\"me:epsilon\"> " + "\n" +
                          "<scalar>" + reactantValues.product1Input.getEpsion() + "</scalar> " + "\n" +
                        "</property> " + "\n" +
                        "<property dictRef=\"me:sigma\"> " + "\n" +
                          "<scalar>" +reactantValues.product1Input.getSigma() + "</scalar>  " + "\n" +
                        "</property> " + "\n" +
                        "<property dictRef=\"me:deltaEDown\"> " + "\n" +
                          "<scalar units=\"cm-1\">" +  reactantValues.product1Input.getDown()+ " </scalar> " + "\n" +
                        "</property> " + "\n" +
                      "</propertyList> " + "\n" +
                      "<me:DOSCMethod>" + reactantValues.product1Input.getDensityOfState()+ "</me:DOSCMethod> " + "\n" +
                    "</molecule> " + "\n";
                 }


            if(reactantValues.getNo_ProductComboSelected().equals("2"))
           {
                  moleculeValue += "   <molecule id=\""+ reactantValues.getProduct_1TextField() + "\">" + "\n" +
                "     <propertyList>"+ "\n" +
                "       <property dictRef=\"me:ZPE\">"+ "\n" +
                "          <scalar units=\""+ measuringUnitValue +"\">" + reactantValues.product2Input.getEnergy_1() + "</scalar>"+ "\n" +
                "           </property>"+ "\n" +
                        "<property dictRef=\"me:rotConsts\"> " + "\n" +
                          "<array units=\"cm-1\">" + reactantValues.product2Input.getRotationalConstant_1()+"</array> " + "\n" +
                        "</property> "+ "\n" +
                        "<property dictRef=\"me:vibFreqs\"> " + "\n" +
                          "<array units=\"cm-1\">" + reactantValues.product2Input.getVibrationalFrequency_1() + "</array> " + "\n" +
                        "</property> "+ "\n" +
                       "<property dictRef=\"me:frequenciesScaleFactor\">" + "\n"+
                        "<scalar>"+   reactantValues.product2Input.getScalingFactor_1() +  "</scalar>" + "\n"+
                       "</property>" + "\n"+
                        "<property dictRef=\"me:symmetryNumber\">" + "\n" +
                          "<scalar>" + reactantValues.product2Input.getSymmetry_1() + "</scalar> "+ "\n" +
                        "</property> " + "\n" +
                        "<property dictRef=\"me:MW\"> " + "\n" +
                          "<scalar units=\"amu\">" + reactantValues.product2Input.getMolecularWeight_1() + "</scalar> " + "\n" +
                        "</property> " + "\n" +
                         "<property dictRef=\"me:spinMultiplicity\"> " + "\n" +
                          "<scalar>" + reactantValues.product2Input.getSpinMultiplicity_1() + "</scalar> " + "\n" +
                        "</property> " + "\n" +
                      "</propertyList> " + "\n" +
                      "<me:DOSCMethod>" + reactantValues.product2Input.getDensityofState_1()+ "</me:DOSCMethod> " + "\n" +
                    "</molecule> " + "\n";

                     // the second recatnat
                    moleculeValue += "   <molecule id=\""+ reactantValues.getProduct_2TextField() + "\">" + "\n" +
                "     <propertyList>"+ "\n" +
                "       <property dictRef=\"me:ZPE\">"+ "\n" +
                "          <scalar units=\""+ measuringUnitValue +"\">" + reactantValues.product2Input.getEnergy_2() + "</scalar>"+ "\n" +
                "           </property>"+ "\n" +
                        "<property dictRef=\"me:rotConsts\"> " + "\n" +
                          "<array units=\"cm-1\">" + reactantValues.product2Input.getRotationalConstant_2()+"</array> " + "\n" +
                        "</property> "+ "\n" +
                        "<property dictRef=\"me:vibFreqs\"> " + "\n" +
                          "<array units=\"cm-1\">" + reactantValues.product2Input.getVibrationalFrequency_2() + "</array> " + "\n" +
                        "</property> "+ "\n" +
                       "<property dictRef=\"me:frequenciesScaleFactor\">" + "\n"+
                        "<scalar>"+   reactantValues.product2Input.getScalingFactor_2() +  "</scalar>" + "\n"+
                       "</property>" + "\n"+
                        "<property dictRef=\"me:symmetryNumber\">" + "\n" +
                          "<scalar>" + reactantValues.product2Input.getSymmetry_2() + "</scalar> "+ "\n" +
                        "</property> " + "\n" +
                        "<property dictRef=\"me:MW\"> " + "\n" +
                          "<scalar units=\"amu\">" + reactantValues.product2Input.getMolecularWeight_2() + "</scalar> " + "\n" +
                        "</property> " + "\n" +
                         "<property dictRef=\"me:spinMultiplicity\"> " + "\n" +
                          "<scalar>" + reactantValues.product2Input.getSpinMultiplicity_2() + "</scalar> " + "\n" +
                        "</property> " + "\n" +
                      "</propertyList> " + "\n" +
                      "<me:DOSCMethod>" + reactantValues.product2Input.getDensityofState_2()+ "</me:DOSCMethod> " + "\n" +
                    "</molecule> " + "\n";
            }

            if(reactantValues.get_RRKM() == true)
            {
                 moleculeValue += "<molecule id=\"" +reactantValues.get_TransitionStateValue() + "\"> " + "\n" +
                          "<propertyList> " + "\n" +
                            "<property dictRef=\"me:ZPE\"> " + "\n" +
                              "<scalar units=\"" +  measuringUnitValue +"\">" + reactantValues.transitionStateInput.getGetEnergy()+"</scalar> " + "\n" +
                            "</property> " + "\n" +
                             "<property dictRef=\"me:rotConsts\"> " + "\n" +
                              "<array units=\"cm-1\">" +reactantValues.transitionStateInput.getRotationalConstants() + "</array> " + "\n" +
                            "</property> " + "\n" +
                             "<property dictRef=\"me:vibFreqs\"> " + "\n" +
                              "<array units=\"cm-1\">"+ reactantValues.transitionStateInput.getVibrationalFrequecy()+"</array> " + "\n" +
                            "</property> " + "\n" +
                            "<property dictRef=\"me:frequenciesScaleFactor\">" + "\n"+
                                "<scalar>"+    reactantValues.transitionStateInput.getScalingFactor() +  "</scalar>" + "\n"+
                           "</property>" + "\n"+
                            "<property dictRef=\"me:symmetryNumber\"> " + "\n" +
                              "<scalar>" +reactantValues.transitionStateInput.getSymmerty()  +"</scalar> " + "\n" +
                            "</property> " + "\n" +
                            "<property dictRef=\"me:MW\"> " + "\n" +
                              "<scalar units=\"amu\">"+ reactantValues.transitionStateInput.getMolecularWeight() +  "</scalar> " + "\n" +
                            "</property> " + "\n" +
                            "<property dictRef=\"me:spinMultiplicity\"> " + "\n" +
                              "<scalar>" + reactantValues.transitionStateInput.getSpinMulticity() +"</scalar> " + "\n" +
                            "</property> " + "\n" +
                            /////////
                            ////    System.out.println("****** " + reactantValues.transitionStateInput.getImaginitiveFrequency() );
                            ///////
                          "</propertyList> " + "\n" +
                          "<me:DOSCMethod>" + reactantValues.transitionStateInput.getDensityOfState()+ "</me:DOSCMethod> " + "\n" +
                        "</molecule> " + "\n";
            }
        }//end of first for loop
        //end of modlule list
         moleculeValue +=       "</moleculeList> " + "\n";

        ReactionPan2 reactantValuesSecond; //this is to go through the reaction list again

        moleculeValue += "<reactionList> " + "\n";  // start of the reaction list
        for(int i =0; i< components.length; i++)
        {
            reactantValuesSecond = (ReactionPan2) components[i];
            moleculeValue+= "<reaction id=\"R"+ (i+1)+ "\"> " + "\n";
            if(reactantValuesSecond.getNo_ReactantComboSelected().equals("1") )
            {
              
                moleculeValue+=   "<reactant> " + "\n" +
                          "<molecule ref=\"" + reactantValuesSecond.getReactant_1TextField()+ "\" me:type=\"" + convertRecProcType(reactantValuesSecond.getReactantCombo1()) + "\" /> " + "\n" +
                          "</reactant> " + "\n";
            }
             if(reactantValuesSecond.getNo_ReactantComboSelected().equals("2") )
            {
                   moleculeValue+=   "<reactant> " + "\n" +
                          "<molecule ref=\"" + reactantValuesSecond.getReactant_1TextField()+ "\" me:type=\"" + convertRecProcType(reactantValuesSecond.getReactantCombo1()) + "\" /> " + "\n" +
                          "</reactant> " + "\n" +
                          "<reactant> " + "\n" +
                          "<molecule ref=\"" + reactantValuesSecond.getReactant_2TextField()+ "\" me:type=\"" + convertRecProcType(reactantValuesSecond.getReactantCombo2()) + "\" /> " + "\n" +
                          "</reactant> " + "\n";
             }

             if(reactantValuesSecond.getNo_ProductComboSelected().equals("1"))
            {
                  moleculeValue+=  "<product> " + "\n" +
                            "<molecule ref=\"" + reactantValuesSecond.getProduct_1TextField()+ "\" me:type=\""+ convertRecProcType(reactantValuesSecond.getProductCombo1()) + "\" /> " + "\n" +
                          "</product> " + "\n";
             }
            
            if(reactantValuesSecond.getNo_ProductComboSelected().equals("2"))
            {
                  moleculeValue+=  "<product> " + "\n" +
                            "<molecule ref=\"" + reactantValuesSecond.getProduct_1TextField()+ "\" me:type=\""+ convertRecProcType(reactantValuesSecond.getProductCombo1()) + "\" /> " + "\n" +
                          "</product> " + "\n" +
                           "<product> " + "\n" +
                            "<molecule ref=\"" + reactantValuesSecond.getProduct_2TextField()+ "\" me:type=\""+ convertRecProcType(reactantValuesSecond.getProductCombo2()) + "\" /> " + "\n" +
                          "</product> " + "\n";
             }

              if(reactantValuesSecond.get_RRKM() == false)
              {
                  moleculeValue+= "<me:MCRCMethod>MesmerILT</me:MCRCMethod> " + "\n" +
                          "<me:preExponential>" +reactantValuesSecond.iltMethodInput.get_ILT_A() + "</me:preExponential> " + "\n" +
                          "<me:excessReactantConc>"  +reactantValuesSecond.iltMethodInput.get_ExcessRegionalConcentration() + "</me:excessReactantConc> " + "\n" +
                          "<me:activationEnergy units=\""+ measuringUnitValue +"\">" + reactantValuesSecond.iltMethodInput.get_ILT_E()+ "</me:activationEnergy> " + "\n" +
                          "<me:nInfinity>" +reactantValuesSecond.iltMethodInput.get_ILT_N() + "</me:nInfinity> " + "\n" +
                        "</reaction> " + "\n";
              }

             if(reactantValuesSecond.get_RRKM() == true)
              {
                  moleculeValue+= "<me:transitionState> " + "\n" +
                            "<molecule ref=\"" + reactantValuesSecond.get_TransitionStateValue() +"\" me:type=\"transitionState\" /> " + "\n" +
                          "</me:transitionState> " + "\n" +
                          "<!--Use zero point energies instead--> " + "\n" +
                          "<me:MCRCMethod>SimpleRRKM</me:MCRCMethod> " + "\n" +
                        "</reaction> " + "\n" ;
             }

          
         
            }// end of second for loop
          moleculeValue+=      "</reactionList> " + "\n";

           moleculeValue+= "<me:conditions> " + "\n" +
                        "<me:bathGas>" + controlObject.get_BufferGas() +"</me:bathGas> " + "\n" +
                        "<me:PTs> " + "\n";
          ArrayList<String> pressure = reactionConditionObject.get_PressureValues();
          ArrayList<String> temperature = reactionConditionObject.get_TemperatureValues();
          for(int i = 0; i < pressure.size() ; i++)
           {
                moleculeValue+= "<me:PTpair me:units=\"" +reactionConditionObject.get_UnitCalculation() + "\" me:P=\"" +pressure.get(i) + "\" me:T=\"" + temperature.get(i)+ "\" /> " + "\n";
//                          "<!--<me:PTpair me:units=\"Torr\" me:P=\"1000\" me:T=\"325.0\" /> " + "\n" +
//                          "<me:PTpair me:units=\"Torr\" me:P=\"1000\" me:T=\"355.0\" />--> " + "\n" +
           }
          moleculeValue+= "</me:PTs> " + "\n" +
                      "</me:conditions> " + "\n";

                      //Told to copy and past by Chem - Mark and Rob but double check again.
          moleculeValue+=  "<me:modelParameters> " + "\n" +
                        "<!--Specify grain size directly...--> " + "\n" +
                        "<me:grainSize units=\"cm-1\">50</me:grainSize> " + "\n" +
                        "<!--...or by the total number of grains " + "\n" +
                            "<me:numberOfGrains> 500 </me:numberOfGrains>--> " + "\n" +
                        "<!--Specify increased energy range " + "\n" +
                            "<me:maxTemperature>6000</me:maxTemperature>--> " + "\n" +
                        "<me:energyAboveTheTopHill>20.</me:energyAboveTheTopHill> " + "\n" +
                      "</me:modelParameters> " + "\n";

          //Get the controlls
          moleculeValue+=   "<me:control> " + "\n";

           if( controlObject.get_calculateRateCoefficientsOnly()  )
            {
                moleculeValue+=  "<me:calculateRateCoefficientsOnly /> " + "\n";
            }
           if(  controlObject.get_printReactionOperatorColumnSums() )
            {
                moleculeValue+=   "<me:printReactionOperatorColumnSums /> " + "\n";
            }
           if(  controlObject.get_printGrainkbE() )
            {
                moleculeValue+= "<me:printGrainkbE /> " + "\n";
            }
          if(     controlObject.get_printReactionOperatorSize() )
            {
                moleculeValue+= "<me:printReactionOperatorSize /> " + "\n";
            }
          if(     controlObject.get_printTunnellingCoefficients()  )
            {
                  moleculeValue+= "<me:printTunnellingCoefficients /> " + "\n";
	   }
         if(     controlObject.get_testRateConstants()  )
	   {
                moleculeValue+= "<me:testRateConstants /> " + "\n";
	   }
         if(     controlObject.get_printCellDOS()  )
	  {
                moleculeValue+= "<me:printCellDOS /> " + "\n";
	  }
         if(     controlObject.get_printGrainBoltzmann()  )
	 {
                moleculeValue+= "<me:printGrainBoltzmann /> " + "\n";
	 }
        if(     controlObject.get_printGrainkfE()  )
	{
            moleculeValue+= "<me:printGrainkfE /> " + "\n";
	}
        if(     controlObject.get_printSpeciesProfile()   )
	{
            moleculeValue+= "<me:printSpeciesProfile /> " + "\n";
	}
        if(     controlObject.get_testDOS()  )
	{
                moleculeValue+= "<me:testDOS /> " + "\n";
	}
        if(     controlObject.get_useTheSameCellNumberForAllConditions() )
	{
            moleculeValue+= "<me:useTheSameCellNumberForAllConditions /> " + "\n";
	}
        if(     controlObject.get_printCellTransitionStateFlux()  )
	{
            moleculeValue+= "<me:printCellTransitionStateFlux /> " + "\n";
	}
        if(     controlObject.get_printGrainDOS()  )
	{
            moleculeValue+= "<me:printGrainDOS /> " + "\n";
	}
        if(     controlObject.get_printGrainTransitionStateFlux()  )
	{
            moleculeValue+= "<me:printGrainTransitionStateFlux /> " + "\n";
	}
        if(     controlObject.get_printTunnelingCoefficients()  )
	{
            moleculeValue+= "<me:printTunnelingCoefficients /> " + "\n";
	}
        if(     controlObject.get_testMicroRates()  )
	{
            moleculeValue+= "<me:testMicroRates /> " + "\n";
	}
                                
        moleculeValue+=  "<me:eigenvalues>" +controlObject.get_Eigenvalues() + "</me:eigenvalues>" + "\n" +
        "</me:control> " + "\n" +
        "</me:mesmer> " + "\n" ;

        xMLStore =""; //set it to empty each time you generate the XML
        xMLStore =  header + "\n"+ moleculeValue;
        System.out.println("XML ******* \n\n\n" + header + "\n"+ moleculeValue);
        System.out.println("XML ******* \n\n\n" + header + "\n"+ moleculeValue);

        saveXMLFile();
    }

    public void saveXMLFile()
    {

        JFileChooser fc = new JFileChooser();
        String outputXMLFIleName = "";

       // fc.setFileSelectionMode(JFileChooser.FILES_AND_DIRECTORIES);
        int returnedValue =  fc.showSaveDialog(mainPanel);
      //  jTextArea1.setText(fc.getSelectedFile().getAbsolutePath());
        if( returnedValue == javax.swing.JFileChooser.APPROVE_OPTION)
        {
            outputXMLFIleName = fc.getSelectedFile().getAbsolutePath();
            System.out.println(" XML file name " + outputXMLFIleName );

            try{
            // Create file
                FileWriter fstream = new FileWriter((outputXMLFIleName));
                BufferedWriter out = new BufferedWriter(fstream);
                out.write(xMLStore);
                //Close the output stream
                out.close();
                }
            catch (Exception e){//Catch exception if any
                System.err.println("Error: " + e.getMessage());
                }
        }


    }
// Commented XML is at the bottom of the page
    
    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        new MESMER_v3().setVisible(true);
    }

}




//
//                     "<molecule id=\"TS_1a\"> " + "\n" +
//                          "<propertyList> " + "\n" +
//                            "<property dictRef=\"me:ZPE\"> " + "\n" +
//                              "<scalar units=\"kJ/mol\">141.4</scalar> " + "\n" +
//                            "</property> " + "\n" +
//                            "<property dictRef=\"me:rotConsts\"> " + "\n" +
//                              "<array units=\"cm-1\">0.33065 0.39313 2.06178</array> " + "\n" +
//                            "</property> " + "\n" +
//                            "<property dictRef=\"me:symmetryNumber\"> " + "\n" +
//                              "<scalar>1.0</scalar> " + "\n" +
//                            "</property> " + "\n" +
//                            "<property dictRef=\"me:vibFreqs\"> " + "\n" +
//                              "<array units=\"cm-1\">332.0 807.0 903.0 1010.0 1098.0 1158.0 1306.0 1407.0 1988.0 3127.0 3167.0</array> " + "\n" +
//                            "</property> " + "\n" +
//                            "<property dictRef=\"me:MW\"> " + "\n" +
//                              "<scalar units=\"amu\">43</scalar> " + "\n" +
//                            "</property> " + "\n" +
//                            "<property dictRef=\"me:spinMultiplicity\"> " + "\n" +
//                              "<scalar>2</scalar> " + "\n" +
//                            "</property> " + "\n" +
//                          "</propertyList> " + "\n" +
//                          "<me:DOSCMethod>ClassicalRotors</me:DOSCMethod> " + "\n" +
//                        "</molecule> " + "\n" +
//                        "<molecule id=\"OCHCH2\"> " + "\n" +
//                          "<propertyList> " + "\n" +
//                            "<property dictRef=\"me:ZPE\"> " + "\n" +
//                              "<scalar units=\"kJ/mol\">-110.5</scalar> " + "\n" +
//                            "</property> " + "\n" +
//                            "<property dictRef=\"me:rotConsts\"> " + "\n" +
//                              "<array units=\"cm-1\">2.247 0.382 0.326</array> " + "\n" +
//                            "</property> " + "\n" +
//                            "<property dictRef=\"me:symmetryNumber\"> " + "\n" +
//                              "<scalar>1.0</scalar> " + "\n" +
//                            "</property> " + "\n" +
//                            "<property dictRef=\"me:vibFreqs\"> " + "\n" +
//                              "<array units=\"cm-1\">443 507 760 975 978 1157 1395 1471 1545 2943 3138 3253</array> " + "\n" +
//                            "</property> " + "\n" +
//                            "<property dictRef=\"me:MW\"> " + "\n" +
//                              "<scalar units=\"amu\">43</scalar> " + "\n" +
//                            "</property> " + "\n" +
//                            "<property dictRef=\"me:spinMultiplicity\"> " + "\n" +
//                              "<scalar>2</scalar> " + "\n" +
//                            "</property> " + "\n" +
//                          "</propertyList> " + "\n" +
//                        "</molecule> " + "\n" +
//                        "<molecule id=\"OH\"> " + "\n" +
//                          "<propertyList> " + "\n" +
//                            "<property dictRef=\"me:ZPE\"> " + "\n" +
//                              "<scalar units=\"kJ/mol\">0.0</scalar> " + "\n" +
//                            "</property> " + "\n" +
//                            "<property dictRef=\"me:rotConsts\"> " + "\n" +
//                              "<array units=\"cm-1\">19.2438</array> " + "\n" +
//                            "</property> " + "\n" +
//                            "<property dictRef=\"me:symmetryNumber\"> " + "\n" +
//                              "<scalar>1</scalar> " + "\n" +
//                            "</property> " + "\n" +
//                            "<property dictRef=\"me:vibFreqs\"> " + "\n" +
//                              "<array units=\"cm-1\">3722.0</array> " + "\n" +
//                            "</property> " + "\n" +
//                            "<property dictRef=\"me:MW\"> " + "\n" +
//                              "<scalar units=\"amu\">17</scalar> " + "\n" +
//                            "</property> " + "\n" +
//                            "<property dictRef=\"me:spinMultiplicity\"> " + "\n" +
//                              "<scalar>2</scalar> " + "\n" +
//                            "</property> " + "\n" +
//                            "<property dictRef=\"me:electronicExcitation\"> " + "\n" +
//                              "<array units=\"cm-1\">139.7</array> " + "\n" +
//                            "</property> " + "\n" +
//                          "</propertyList> " + "\n" +
//                          "<me:DOSCMethod>ClassicalRotors</me:DOSCMethod> " + "\n" +
//                        "</molecule> " + "\n" +
//                        "<molecule id=\"acetylene\"> " + "\n" +
//                          "<propertyList> " + "\n" +
//                            "<property dictRef=\"me:ZPE\"> " + "\n" +
//                              //"<scalar units=\"kJ/mol\">130.7</scalar> \"
//                             "<scalar units=\"kJ/mol\">130.7</scalar> " + "\n" +
//                            "</property> " + "\n" +
//                            "<property dictRef=\"me:rotConsts\"> " + "\n" +
//                              "<array units=\"cm-1\">1.19205</array> " + "\n" +
//                            "</property> " + "\n" +
//                            "<property dictRef=\"me:symmetryNumber\"> " + "\n" +
//                              "<scalar>2</scalar> " + "\n" +
//                            "</property> " + "\n" +
//                            "<property dictRef=\"me:vibFreqs\"> " + "\n" +
//                              "<array units=\"cm-1\">664.0 664.0 766.0 766.0 2066.0 3408.0 3509.0</array> " + "\n" +
//                            "</property> " + "\n" +
//                            "<property dictRef=\"me:MW\"> " + "\n" +
//                              "<scalar units=\"amu\">26</scalar> " + "\n" +
//                            "</property> " + "\n" +
//                            "<property dictRef=\"me:spinMultiplicity\"> " + "\n" +
//                              "<scalar>1</scalar> " + "\n" +
//                            "</property> " + "\n" +
//                          "</propertyList> " + "\n" +
//                          "<me:DOSCMethod>ClassicalRotors</me:DOSCMethod> " + "\n" +
//                        "</molecule> " + "\n" +
//                      "</moleculeList> " + "\n" +
//                      "<reactionList> " + "\n" +
//
//
//                      "<reaction id=\"R1\"> " + "\n" +
//                          "<reactant> " + "\n" +
//                            "<molecule ref=\"OH\" me:type=\"deficientReactant\" /> " + "\n" +
//                          "</reactant> " + "\n" +
//                          "<reactant> " + "\n" +
//                            "<molecule ref=\"acetylene\" me:type=\"excessReactant\" /> " + "\n" +
//                          "</reactant> " + "\n" +
//                          "<product> " + "\n" +
//                            "<molecule ref=\"OH_acetylene_adduct\" me:type=\"modelled\" /> " + "\n" +
//                          "</product> " + "\n" +
//
//
//
//
//                          "<me:MCRCMethod>MesmerILT</me:MCRCMethod> " + "\n" +
//                          "<me:preExponential>7.00e-12</me:preExponential> " + "\n" +
//                          "<me:excessReactantConc>1.00e17</me:excessReactantConc> " + "\n" +
//                          "<me:activationEnergy units=\"kJ/mol\">5.375</me:activationEnergy> " + "\n" +
//                          "<me:nInfinity>0.0</me:nInfinity> " + "\n" +
//                        "</reaction> " + "\n" +
//                        "<reaction id=\"R2\"> " + "\n" +
//                          "<reactant> " + "\n" +
//                            "<molecule ref=\"OH_acetylene_adduct\" me:type=\"reactant\" /> " + "\n" +
//                          "</reactant> " + "\n" +
//                          "<product> " + "\n" +
//                            "<molecule ref=\"OCHCH2\" me:type=\"sink\" /> " + "\n" +
//                          "</product> " + "\n" +
//                          "<me:transitionState> " + "\n" +
//                            "<molecule ref=\"TS_1a\" me:type=\"transitionState\" /> " + "\n" +
//                          "</me:transitionState> " + "\n" +
//                          "<!--Use zero point energies instead--> " + "\n" +
//                          "<me:MCRCMethod>SimpleRRKM</me:MCRCMethod> " + "\n" +
//                        "</reaction> " + "\n" +
//
//
//                      "</reactionList> " + "\n" +
//                      "<me:conditions> " + "\n" +
//                        "<me:bathGas>N2</me:bathGas> " + "\n" +
//                        "<me:PTs> " + "\n" +
//                          "<me:PTpair me:units=\"Torr\" me:P=\"1000\" me:T=\"298.0\" /> " + "\n" +
//                          "<!--<me:PTpair me:units=\"Torr\" me:P=\"1000\" me:T=\"325.0\" /> " + "\n" +
//                          "<me:PTpair me:units=\"Torr\" me:P=\"1000\" me:T=\"355.0\" />--> " + "\n" +
//                        "</me:PTs> " + "\n" +
//                      "</me:conditions> " + "\n" +
//                      "<me:modelParameters> " + "\n" +
//                        "<!--Specify grain size directly...--> " + "\n" +
//                        "<me:grainSize units=\"cm-1\">50</me:grainSize> " + "\n" +
//                        "<!--...or by the total number of grains " + "\n" +
//                            "<me:numberOfGrains> 500 </me:numberOfGrains>--> " + "\n" +
//                        "<!--Specify increased energy range " + "\n" +
//                            "<me:maxTemperature>6000</me:maxTemperature>--> " + "\n" +
//                        "<me:energyAboveTheTopHill>20.</me:energyAboveTheTopHill> " + "\n" +
//                      "</me:modelParameters> " + "\n" +
//
//
//
//                      "<me:control> " + "\n" +
//                        "<me:testDOS /> " + "\n" +
//                        "<me:testMicroRates /> " + "\n" +
//                        "<!--<me:printGrainDOS />--> " + "\n" +
//                        "<!--<me:printCellDOS />--> " + "\n" +
//                        "<!--<me:printReactionOperatorColumnSums />--> " + "\n" +
//                        "<me:printGrainkfE /> " + "\n" +
//                        "<!--<me:printGrainBoltzmann />--> " + "\n" +
//                        "<me:printGrainkbE /> " + "\n" +
//                        "<me:printSpeciesProfile /> " + "\n" +
//                        "<me:testRateConstants /> " + "\n" +
//                        "<me:eigenvalues>0</me:eigenvalues> " + "\n" +
//                      "</me:control> " + "\n" +
//                    "</me:mesmer> " + "\n" ;



   /////////////////////////////



//            }
//      }
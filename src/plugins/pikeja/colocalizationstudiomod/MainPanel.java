package plugins.pikeja.colocalizationstudiomod;

import javax.swing.JPanel;
import javax.swing.JButton;
import javax.swing.JLabel;
import javax.swing.border.TitledBorder;
import javax.swing.SwingConstants;
import java.awt.GridBagLayout;
import java.awt.GridBagConstraints;
import java.awt.Insets;
import javax.swing.JComboBox;

public class MainPanel extends JPanel {
	
	public JButton buttonStartComputation;
	public JButton roiSetAButton;
	public JLabel labelSetA;
	public JButton roiSetBButton;
	public JLabel labelSetB;
	public JComboBox comboBoxColocMethod;

	/**
	 * Create the panel.
	 */
	public MainPanel() {
		GridBagLayout gridBagLayout = new GridBagLayout();
		gridBagLayout.columnWidths = new int[]{0, 0};
		gridBagLayout.rowHeights = new int[]{0, 0, 0, 0, 0};
		gridBagLayout.columnWeights = new double[]{1.0, Double.MIN_VALUE};
		gridBagLayout.rowWeights = new double[]{1.0, 1.0, 1.0, 1.0, Double.MIN_VALUE};
		setLayout(gridBagLayout);
		
		JPanel panel_1 = new JPanel();
		panel_1.setBorder(new TitledBorder(null, "ROI Set A", TitledBorder.LEADING, TitledBorder.TOP, null, null));
		GridBagConstraints gbc_panel_1 = new GridBagConstraints();
		gbc_panel_1.anchor = GridBagConstraints.NORTH;
		gbc_panel_1.fill = GridBagConstraints.HORIZONTAL;
		gbc_panel_1.insets = new Insets(0, 0, 5, 0);
		gbc_panel_1.gridx = 0;
		gbc_panel_1.gridy = 0;
		add(panel_1, gbc_panel_1);
		GridBagLayout gbl_panel_1 = new GridBagLayout();
		gbl_panel_1.columnWidths = new int[]{0, 0};
		gbl_panel_1.rowHeights = new int[]{23, 14, 0};
		gbl_panel_1.columnWeights = new double[]{1.0, Double.MIN_VALUE};
		gbl_panel_1.rowWeights = new double[]{0.0, 0.0, Double.MIN_VALUE};
		panel_1.setLayout(gbl_panel_1);
		
		roiSetAButton = new JButton("Affect current selected ROI(s) to set A");
		GridBagConstraints gbc_roiSetAButton = new GridBagConstraints();
		gbc_roiSetAButton.fill = GridBagConstraints.HORIZONTAL;
		gbc_roiSetAButton.insets = new Insets(0, 0, 5, 0);
		gbc_roiSetAButton.gridx = 0;
		gbc_roiSetAButton.gridy = 0;
		panel_1.add(roiSetAButton, gbc_roiSetAButton);
		
		labelSetA = new JLabel("**");
		labelSetA.setHorizontalAlignment(SwingConstants.CENTER);
		GridBagConstraints gbc_labelSetA = new GridBagConstraints();
		gbc_labelSetA.fill = GridBagConstraints.HORIZONTAL;
		gbc_labelSetA.gridx = 0;
		gbc_labelSetA.gridy = 1;
		panel_1.add(labelSetA, gbc_labelSetA);
		
		JPanel panel = new JPanel();
		panel.setBorder(new TitledBorder(null, "ROI Set B", TitledBorder.LEADING, TitledBorder.TOP, null, null));
		GridBagConstraints gbc_panel = new GridBagConstraints();
		gbc_panel.fill = GridBagConstraints.BOTH;
		gbc_panel.insets = new Insets(0, 0, 5, 0);
		gbc_panel.gridx = 0;
		gbc_panel.gridy = 1;
		add(panel, gbc_panel);
		GridBagLayout gbl_panel = new GridBagLayout();
		gbl_panel.columnWidths = new int[]{0, 0};
		gbl_panel.rowHeights = new int[]{23, 14, 0};
		gbl_panel.columnWeights = new double[]{1.0, Double.MIN_VALUE};
		gbl_panel.rowWeights = new double[]{0.0, 0.0, Double.MIN_VALUE};
		panel.setLayout(gbl_panel);
		
		roiSetBButton = new JButton("Affect current slected ROI(s) to set B");
		GridBagConstraints gbc_roiSetBButton = new GridBagConstraints();
		gbc_roiSetBButton.fill = GridBagConstraints.HORIZONTAL;
		gbc_roiSetBButton.insets = new Insets(0, 0, 5, 0);
		gbc_roiSetBButton.gridx = 0;
		gbc_roiSetBButton.gridy = 0;
		panel.add(roiSetBButton, gbc_roiSetBButton);
		
		labelSetB = new JLabel("**");
		labelSetB.setHorizontalAlignment(SwingConstants.CENTER);
		GridBagConstraints gbc_labelSetB = new GridBagConstraints();
		gbc_labelSetB.fill = GridBagConstraints.HORIZONTAL;
		gbc_labelSetB.gridx = 0;
		gbc_labelSetB.gridy = 1;
		panel.add(labelSetB, gbc_labelSetB);
		
		JPanel panel_3 = new JPanel();
		panel_3.setBorder(new TitledBorder(null, "Colocalization options", TitledBorder.LEADING, TitledBorder.TOP, null, null));
		GridBagConstraints gbc_panel_3 = new GridBagConstraints();
		gbc_panel_3.fill = GridBagConstraints.BOTH;
		gbc_panel_3.insets = new Insets(0, 0, 5, 0);
		gbc_panel_3.gridx = 0;
		gbc_panel_3.gridy = 2;
		add(panel_3, gbc_panel_3);
		GridBagLayout gbl_panel_3 = new GridBagLayout();
		gbl_panel_3.columnWidths = new int[]{0, 0};
		gbl_panel_3.rowHeights = new int[]{0, 0};
		gbl_panel_3.columnWeights = new double[]{1.0, Double.MIN_VALUE};
		gbl_panel_3.rowWeights = new double[]{0.0, Double.MIN_VALUE};
		panel_3.setLayout(gbl_panel_3);
		
		comboBoxColocMethod = new JComboBox();
		GridBagConstraints gbc_comboBoxColocMethod = new GridBagConstraints();
		gbc_comboBoxColocMethod.fill = GridBagConstraints.HORIZONTAL;
		gbc_comboBoxColocMethod.gridx = 0;
		gbc_comboBoxColocMethod.gridy = 0;
		panel_3.add(comboBoxColocMethod, gbc_comboBoxColocMethod);
		
		JPanel panel_2 = new JPanel();
		panel_2.setBorder(new TitledBorder(null, "Computation", TitledBorder.LEADING, TitledBorder.TOP, null, null));
		GridBagConstraints gbc_panel_2 = new GridBagConstraints();
		gbc_panel_2.anchor = GridBagConstraints.NORTH;
		gbc_panel_2.fill = GridBagConstraints.HORIZONTAL;
		gbc_panel_2.gridx = 0;
		gbc_panel_2.gridy = 3;
		add(panel_2, gbc_panel_2);
		GridBagLayout gbl_panel_2 = new GridBagLayout();
		gbl_panel_2.columnWidths = new int[]{0, 0};
		gbl_panel_2.rowHeights = new int[]{23, 0};
		gbl_panel_2.columnWeights = new double[]{1.0, Double.MIN_VALUE};
		gbl_panel_2.rowWeights = new double[]{0.0, Double.MIN_VALUE};
		panel_2.setLayout(gbl_panel_2);
		
		buttonStartComputation = new JButton("Start computation");
		GridBagConstraints gbc_buttonStartComputation = new GridBagConstraints();
		gbc_buttonStartComputation.fill = GridBagConstraints.HORIZONTAL;
		gbc_buttonStartComputation.gridx = 0;
		gbc_buttonStartComputation.gridy = 0;
		panel_2.add(buttonStartComputation, gbc_buttonStartComputation);

	}

}

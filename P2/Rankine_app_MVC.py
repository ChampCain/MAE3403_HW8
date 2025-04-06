#region imports
import sys
from PyQt5 import QtWidgets as qtw
from PyQt5 import QtCore as qtc
from Rankine_GUI import Ui_Form
from Rankine_Classes_MVC import rankineController
from UnitConversions import UnitConverter as UC
#these imports are necessary for drawing a matplot lib graph on my GUI
#no simple widget for this exists in QT Designer, so I have to add the widget in code.
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg, NavigationToolbar2QT
from matplotlib.figure import Figure
#endregion

#region class definitions
class MainWindow(qtw.QWidget, Ui_Form):
    def __init__(self):
        """
        MainWindow constructor
        """
        super().__init__()
        # Main UI setup
        self.setupUi(self)
        self.AssignSlots()
        self.MakeCanvas()

        # Create lists of input and display widgets
        self.input_widgets = [self.rb_SI, self.le_PHigh, self.le_PLow,
                              self.le_TurbineInletCondition, self.rdo_Quality,
                              self.le_TurbineEff, self.cmb_XAxis, self.cmb_YAxis,
                              self.chk_logX, self.chk_logY]

        self.display_widgets = [self.lbl_PHigh, self.lbl_PLow, self.lbl_SatPropLow,
                                self.lbl_SatPropHigh, self.lbl_TurbineInletCondition,
                                self.lbl_H1, self.lbl_H1Units, self.lbl_H2,
                                self.lbl_H2Units, self.lbl_H3, self.lbl_H3Units,
                                self.lbl_H4, self.lbl_H4Units, self.lbl_TurbineWork,
                                self.lbl_TurbineWorkUnits, self.lbl_PumpWork,
                                self.lbl_PumpWorkUnits, self.lbl_HeatAdded,
                                self.lbl_HeatAddedUnits, self.lbl_ThermalEfficiency,
                                self.canvas, self.figure, self.ax]

        # Instantiate controller
        self.RC = rankineController(self.input_widgets, self.display_widgets)

        # Initialize pressure displays and turbine inlet
        self.setNewPHigh()
        self.setNewPLow()
        self.SelectQualityOrTHigh()  # Sets proper initial value

        self.Calculate()  # Perform initial calculation

        # Store coordinates from last position on graph
        self.oldXData = 0.0
        self.oldYData = 0.0

        self.show()

    def AssignSlots(self):
        """
        Setup signals and slots for my program
        """
        self.btn_Calculate.clicked.connect(self.Calculate)
        self.rdo_Quality.clicked.connect(self.SelectQualityOrTHigh)
        self.rdo_THigh.clicked.connect(self.SelectQualityOrTHigh)
        self.le_PHigh.editingFinished.connect(self.setNewPHigh)
        self.le_PLow.editingFinished.connect(self.setNewPLow)
        self.rb_SI.clicked.connect(self.SetUnits)
        self.rb_English.clicked.connect(self.SetUnits)
        self.cmb_XAxis.currentIndexChanged.connect(self.SetPlotVariables)
        self.cmb_YAxis.currentIndexChanged.connect(self.SetPlotVariables)
        self.chk_logX.toggled.connect(self.SetPlotVariables)
        self.chk_logY.toggled.connect(self.SetPlotVariables)

    def MakeCanvas(self):
        """
        Create a place to make graph of Rankine cycle
        """
        self.figure = Figure(figsize=(1, 1), tight_layout=True, frameon=True)
        self.canvas = FigureCanvasQTAgg(self.figure)
        self.ax = self.figure.add_subplot()
        self.Layout_Plot.addWidget(NavigationToolbar2QT(self.canvas, self))
        self.Layout_Plot.addWidget(self.canvas)
        self.canvas.mpl_connect("motion_notify_event", self.mouseMoveEvent_Canvas)

    def mouseMoveEvent_Canvas(self, event):
        self.oldXData = event.xdata if event.xdata is not None else self.oldXData
        self.oldYData = event.ydata if event.ydata is not None else self.oldYData
        sUnit = 'kJ/(kg*K)' if self.rb_SI.isChecked() else 'BTU/(lb*R)'
        TUnit = 'C' if self.rb_SI.isChecked() else 'F'
        self.setWindowTitle('s:{:0.2f} {}, T:{:0.2f} {}'.format(self.oldXData, sUnit, self.oldYData, TUnit))

    def Calculate(self):
        self.RC.updateModel()

    def SelectQualityOrTHigh(self):
        self.RC.selectQualityOrTHigh()

    def SetPlotVariables(self):
        self.RC.updatePlot()

    def SetUnits(self):
        self.RC.updateUnits()

    def setNewPHigh(self):
        """
        Update the saturated properties display when P High changes
        """
        try:
            SI = self.rb_SI.isChecked()
            # Convert input to kPa (steam tables use kPa)
            PCF = 100 if SI else UC.psi_to_kpa  # bar->kPa or psi->kPa
            PHigh_text = self.le_PHigh.text()

            if not PHigh_text:
                return

            PHigh = float(PHigh_text) * PCF

            # Validate pressure range (typical steam tables go up to ~220 bar)
            if PHigh > 22000:  # 220 bar = 22000 kPa
                print(f"Warning: Pressure {PHigh} kPa is too high for steam tables")
                return

            # Get saturation properties
            satPHigh = self.RC.Model.steam.getsatProps_p(PHigh)
            if satPHigh is None:
                print(f"Could not get saturation properties for P={PHigh} kPa")
                return

            # Update saturation properties display
            self.lbl_SatPropHigh.setText(satPHigh.getTextOutput(SI=SI))

            # If in T High mode, update the turbine inlet temperature
            if self.rdo_THigh.isChecked():
                Tsat = satPHigh.tsat
                if not SI:
                    Tsat = UC.C_to_F(Tsat)
                self.le_TurbineInletCondition.setText(f"{Tsat:.2f}")
        except ValueError:
            print("Invalid pressure value entered")
        except Exception as e:
            print(f"Error in setNewPHigh: {e}")

    def setNewPLow(self):
        """
        Update the saturated properties display when P Low changes
        """
        try:
            SI = self.rb_SI.isChecked()
            # Convert input to kPa (steam tables use kPa)
            PCF = 100 if SI else UC.psi_to_kpa  # bar->kPa or psi->kPa
            PLow_text = self.le_PLow.text()

            if not PLow_text:
                return

            PLow = float(PLow_text) * PCF

            # Get saturation properties
            satPLow = self.RC.Model.steam.getsatProps_p(PLow)
            if satPLow is None:
                print(f"Could not get saturation properties for P={PLow} kPa")
                return

            # Update saturation properties display
            self.lbl_SatPropLow.setText(satPLow.getTextOutput(SI=SI))
        except ValueError:
            print("Invalid pressure value entered")
        except Exception as e:
            print(f"Error in setNewPLow: {e}")

#endregion

#if this module is being imported, this won't run. If it is the main module, it will run.
if __name__== '__main__':
    app = qtw.QApplication(sys.argv)
    mw = MainWindow()
    mw.setWindowTitle('Rankine calculator')
    sys.exit(app.exec())
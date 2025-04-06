#region imorts
import numpy as np
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from matplotlib.figure import Figure
import PyQt5.QtWidgets as qtw

# importing from previous work on least squares fit
from LeastSquares import LeastSquaresFit_Class
#endregion

#region class definitions
class Pump_Model():
    """
    This is the pump model.  It just stores data.
    """
    def __init__(self): #pump class constructor
        #create some class variables for storing information
        self.PumpName = ""
        self.FlowUnits = ""
        self.HeadUnits = ""

        # place to store data from file
        self.FlowData = np.array([])
        self.HeadData = np.array([])
        self.EffData = np.array([])

        # place to store coefficients for cubic fits
        self.HeadCoefficients = np.array([])
        self.EfficiencyCoefficients = np.array([])

        # create two instances (objects) of least squares class
        self.LSFitHead=LeastSquaresFit_Class()
        self.LSFitEff=LeastSquaresFit_Class()

class Pump_Controller():
    def __init__(self):
        self.Model = Pump_Model()
        self.View = Pump_View()
    
    #region functions to modify data of the model
    def ImportFromFile(self, data):
        """
        This processes the list of strings in data to build the pump model
        :param data: List of strings containing pump data
        :return: None
        """
        self.Model.PumpName = data[0].strip()  # First line is pump name
        # data[1] is the units line
        L = data[2].split()
        self.Model.FlowUnits = L[0]
        self.Model.HeadUnits = L[1]

        # extracts flow, head and efficiency data and calculates coefficients
        self.SetData(data[3:])
        self.updateView()

    def SetData(self, data):
        '''
        Expects three columns of data in an array of strings with space delimiter
        Parse line and build arrays.
        :param data: List of strings containing numerical data
        :return: None
        '''
        # erase existing data
        self.Model.FlowData = np.array([])
        self.Model.HeadData = np.array([])
        self.Model.EffData = np.array([])

        # parse new data
        for L in data:
            Cells = L.split()  # split the line into an array of strings
            self.Model.FlowData = np.append(self.Model.FlowData, float(Cells[0].strip()))  # first column is flow
            self.Model.HeadData = np.append(self.Model.HeadData, float(Cells[1].strip()))  # second column is head
            self.Model.EffData = np.append(self.Model.EffData, float(Cells[2].strip()))  # third column is efficiency

        # call least square fit for head and efficiency
        self.LSFit()
        
    def LSFit(self):
        '''Fit cubic polynomial using Least Squares'''
        self.Model.LSFitHead.x=self.Model.FlowData
        self.Model.LSFitHead.y=self.Model.HeadData
        self.Model.LSFitHead.LeastSquares(3) #calls LeastSquares function of LSFitHead object

        self.Model.LSFitEff.x=self.Model.FlowData
        self.Model.LSFitEff.y=self.Model.EffData
        self.Model.LSFitEff.LeastSquares(3) #calls LeastSquares function of LSFitEff object
    #endregion

    #region functions interacting with view
    def setViewWidgets(self, w):
        self.View.setViewWidgets(w)

    def updateView(self):
        self.View.updateView(self.Model)
    #endregion
class Pump_View():
    def __init__(self):
        """
        In this constructor, I create some QWidgets as placeholders until they get defined later.
        """
        self.LE_PumpName=qtw.QLineEdit()
        self.LE_FlowUnits=qtw.QLineEdit()
        self.LE_HeadUnits=qtw.QLineEdit()
        self.LE_HeadCoefs=qtw.QLineEdit()
        self.LE_EffCoefs=qtw.QLineEdit()
        self.ax=None
        self.canvas=None

    def updateView(self, Model):
        """
        Put model parameters in the widgets.
        :param Model:
        :return:
        """
        self.LE_PumpName.setText(Model.PumpName)
        self.LE_FlowUnits.setText(Model.FlowUnits)
        self.LE_HeadUnits.setText(Model.HeadUnits)
        self.LE_HeadCoefs.setText(Model.LSFitHead.GetCoeffsString())
        self.LE_EffCoefs.setText(Model.LSFitEff.GetCoeffsString())
        self.DoPlot(Model)

    def DoPlot(self, Model):
        """
        Create a black and white plot with customized styling and legend placement
        """
        # Get fit data - quadratic for head (power=2), cubic for efficiency (power=3)
        headx, heady, headRSq = Model.LSFitHead.GetPlotInfo(2, npoints=500)
        effx, effy, effRSq = Model.LSFitEff.GetPlotInfo(3, npoints=500)

        axes = self.ax
        axes.clear()  # clear any previous plot

        # Plot Head Data (left y-axis) - Hollow black circles
        head_data_line = axes.plot(Model.FlowData, Model.HeadData, 'ko',
                                   markersize=6, markerfacecolor='none',
                                   markeredgewidth=1.5, label='Head')[0]

        # Plot Head Fit - Dashed black line
        head_fit_line = axes.plot(headx, heady, 'k--', linewidth=1.5,
                                  label=f'Head (R²={headRSq:.3f})')[0]

        axes.set_xlabel(f'Flow Rate ({Model.FlowUnits})', fontsize=10)
        axes.set_ylabel(f'Head ({Model.HeadUnits})', color='k', fontsize=10)
        axes.tick_params(axis='y', labelcolor='k')

        # Create second y-axis for Efficiency
        ax2 = axes.twinx()

        # Plot Efficiency Data - Hollow black triangles
        eff_data_line = ax2.plot(Model.FlowData, Model.EffData, 'k^',
                                 markersize=6, markerfacecolor='none',
                                 markeredgewidth=1.5, label='Efficiency')[0]

        # Plot Efficiency Fit - Dotted black line
        eff_fit_line = ax2.plot(effx, effy, 'k:', linewidth=2,
                                label=f'Efficiency (R²={effRSq:.3f})')[0]

        ax2.set_ylabel('Efficiency (%)', color='k', fontsize=10)
        ax2.tick_params(axis='y', labelcolor='k')

        # Set title
        axes.set_title(f'Pump Curve: {Model.PumpName}', fontsize=12)

        # Create separate legends
        # Left legend for Head (upper left)
        head_legend = axes.legend(handles=[head_data_line, head_fit_line],
                                  loc='upper left', fontsize=8,
                                  framealpha=1, edgecolor='k')

        # Add the legend manually to the axes
        axes.add_artist(head_legend)

        # Right legend for Efficiency (upper right)
        eff_legend = ax2.legend(handles=[eff_data_line, eff_fit_line],
                                loc='upper right', fontsize=8,
                                framealpha=1, edgecolor='k')

        # Add grid for better readability
        axes.grid(True, linestyle='--', alpha=0.7, color='k')

        # Auto-scale axes (removed hard boundaries)
        axes.relim()
        axes.autoscale_view()
        ax2.relim()
        ax2.autoscale_view()

        self.canvas.draw()

    def setViewWidgets(self, w):
        self.LE_PumpName, self.LE_FlowUnits, self.LE_HeadUnits, self.LE_HeadCoefs, self.LE_EffCoefs, self.ax, self.canvas = w
#endregion


from psychopy import prefs
prefs.hardware['audioLib'] =  ['PTB', 'sounddevice', 'pyo', 'pygame']
from psychopy import core, visual, gui, data, event, sound
from psychopy.tools.filetools import fromFile, toFile
import numpy as np
import numpy, random
from scipy import stats


#degrees of visual angle:
#please calibrate PsychoPy to match the monitor which it is using through Tools-->Monitor Center

#create window & shapes
mywin=visual.Window(size=[800,600], monitor="testMonitor", units="deg") 
square = visual.Rect(win=mywin, units ='deg', size=[3,3], fillColor='white') #tried to adjust for degree of visual angle here and right below
gabor = visual.GratingStim(mywin, units ='deg', size=[3,3], tex='sin', mask='gauss', sf=5, name='gabor', autoLog=False)

beep = sound.Sound(value='C', secs=0.5, octave=4, stereo=- 1, volume=1.0, loops=0, sampleRate=None, blockSize=128, preBuffer=- 1, hamming=True, startTime=0, stopTime=- 1, name='', autoLog=True)
errorMsg = visual.TextStim(mywin, text='error', color='red')

trialClock = core.Clock()
timer = core.CountdownTimer(.5)

pResponses = []
# create the experiment/trial? handler
experiment = data.ExperimentHandler(dataFileName="results")
trials = data.TrialHandler(pResponses, 100, method='random')
trials.data.addDataType('stimType')
trials.data.addDataType('stimSide')
trials.data.addDataType('response')
trials.data.addDataType('responseRT')
trials.data.addDataType('error')


# display instructions and wait
message1 = visual.TextStim(mywin, pos=[0,+3],text='Hit the space key when ready.')
message2 = visual.TextStim(mywin, pos=[0,-3],
    text="If the square displayed is a Gabor, press the 'z' key when it is on the right, and the 'm' key if it is on the left. However, if this square is white, press the key that is on the same side as the white square.")
message1.draw()
message2.draw()
mywin.flip()#to show our newly drawn 'stimuli'
#pause until there's a keypress
event.waitKeys(keyList='space')

orderList=[i for i in range(100)]
for trial in trials:  # will continue the trials until it terminates!
        # set location of stimuli
        trialClock.reset()
        timer.reset()
        error = 0
        rand=np.random.choice(orderList)
        if rand%2==0:
            shape=gabor
            shapeType=1
            oriList=[i for i in range(359)]
            shape.ori=np.random.choice(oriList)
            orderList.remove(rand)
        else:
            shape=square
            shapeType=0
            orderList.remove(rand)
        #shape=np.random.choice([gabor, square], p=[.5,.5])
        shapeSide= np.random.choice([-1,1])  # will be either +1(right) or -1(left)
        shape.setPos([8*shapeSide, 0])
        while timer.getTime()>0:
            shape.draw()
            mywin.flip()
        mywin.flip(clearBuffer=True)
        # get response
        thisResp=None
        while thisResp==None:
            allKeys=event.getKeys()
            for thisKey in allKeys:
                if thisKey=='m':
                    thisResp = 1 
                    if (shapeType==0 and shapeSide==-1) or (shapeType==1 and shapeSide==1):
                        beep.play()
                        errorMsg.setPos([8*shapeSide, 0])
                        shape.draw()
                        errorMsg.draw()
                        mywin.flip()
                        error = 1
                elif thisKey=='z':
                    thisResp = -1
                    if (shapeType==0 and shapeSide==1) or (shapeType==1 and shapeSide==-1):
                        beep.play()
                        errorMsg.setPos([8*shapeSide, 0])
                        shape.draw()
                        errorMsg.draw()
                        mywin.flip()
                        error = 1
                elif thisKey in ['q', 'escape']:
                    core.quit()  # abort experiment
        core.wait(.3)
        event.clearEvents()  # clear other (eg mouse) events - they clog the buffer
        mywin.flip(clearBuffer=True)
        # add the data to the trial to go to the next iteration
        pResponses.append({'stimType':shapeType,'stimSide':shapeSide, 'response':thisResp, 'responseRT':trialClock.getTime(), 'error':error})
        trials.data.add('stimType',shapeType)
        trials.data.add('stimSide', shapeSide)
        trials.data.add('response', thisResp)
        trials.data.add('responseRT',trialClock.getTime())
        trials.data.add('error',error)
        
        experiment.addData('stimType', shapeType)
        experiment.addData('stimSide', shapeSide)
        experiment.addData('response', thisResp)
        experiment.addData('responseRT',trialClock.getTime())
        experiment.addData('error',error)
        experiment.nextEntry()
        core.wait(1.5)

print (pResponses)
experiment.close()
#dataFile.close()
#staircase.saveAsPickle(fileName) 

#win.close()
#core.quit()
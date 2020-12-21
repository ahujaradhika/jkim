from psychopy import core, visual, gui, data, event, sound
from psychopy.tools.filetools import fromFile, toFile
import numpy, random
from scipy import stats

#degrees of visual angle:
#please calibrate PsychoPy to match the monitor which it is using through Tools-->Monitor Center

#create window & shapes
mywin=visual.Window(size=[800,600], monitor="testMonitor", units="deg")
shape = visual.Rect(win=mywin, units ='deg', size=[3,3], fillColor='white')
#gabor = visual.GratingStim(win, tex='sin', mask='gauss', sf=5, name='gabor', autoLog=False)

beep = sound.Sound(value='C', secs=0.5, octave=4, stereo=- 1, volume=1.0, loops=0, sampleRate=None, blockSize=128, preBuffer=- 1, hamming=True, startTime=0, stopTime=- 1, name='', autoLog=True)
errorMsg = visual.TextStim(mywin, text='error', color='red')

trialClock = core.Clock()

pResponses = []
# create the experiment/trial? handler
trials = data.TrialHandler(pResponses, 3, method='random')
trials.data.addDataType('stimID')
trials.data.addDataType('response')
trials.data.addDataType('responseRT')

# display instructions and wait
message1 = visual.TextStim(mywin, pos=[0,+3],text='Hit the space key when ready.')
message2 = visual.TextStim(mywin, pos=[0,-3],
    text="If the square displayed is a Gabor, press the 'z' key when it is on the right, and the 'm' key if it is on the left. However, if this square is white, press the key that is on the same side as the white square.")
message1.draw()
message2.draw()
mywin.flip()#to show our newly drawn 'stimuli'
#pause until there's a keypress
event.waitKeys(keyList='space')


for trial in trials:  # will continue the trials until it terminates!
        # set location of stimuli
        trialClock.reset()
        shapeSide= random.choice([-1,1])  # will be either +1(right) or -1(left)
        shape.setPos([8*shapeSide, 0])
        # if trialClock.getTime()<=1.0: #<---- THIS LINE NOT WORKING HERE
        shape.draw()
        mywin.flip()
        #elif trialClock.getTime()>1.0: <---- THIS LINE NOT WORKING HERE
            #mywin.flip(clearBuffer=True)
            
        # get response
        thisResp=None
        while thisResp==None:
            allKeys=event.waitKeys()
            for thisKey in allKeys:
                if thisKey=='m':
                    beep.play()
                    thisResp = 1 
                    if shapeSide==-1:
                        beep.play() #<----- THIS ONLY PLAYS IN FIRST TRIAL, NO SOUNDS PLAY AFTER FIRST TRIAL
                        errorMsg.setPos([8*shapeSide, 0])
                        #if trialClock.getTime()<=1.0: <---- THIS LINE NOT WORKING HERE 
                        shape.draw()
                        errorMsg.draw()
                        mywin.flip()
                        core.wait(.3)
                elif thisKey=='z':
                    thisResp = -1
                    if shapeSide==1:
                        beep.play() #<----- THIS ONLY PLAYS IN FIRST TRIAL, NO SOUNDS PLAY AFTER FIRST TRIAL
                        errorMsg.setPos([8*shapeSide, 0])
                        #if trialClock.getTime()<=1.0: <---- THIS LINE NOT WORKING HERE
                        shape.draw()
                        errorMsg.draw()
                        mywin.flip()
                        core.wait(.3)
                elif thisKey in ['q', 'escape']:
                    core.quit()  # abort experiment
        event.clearEvents()  # clear other (eg mouse) events - they clog the buffer
        mywin.flip(clearBuffer=True)
        # add the data to the trial to go to the next iteration
        pResponses.append({'stimID':shapeSide, 'response':thisResp, 'responseRT':trialClock.getTime()})
        trials.data.add('stimID', shapeSide)
        trials.data.add('response', thisResp)
        trials.data.add('responseRT',trialClock.getTime())
        core.wait(1.5)

print (pResponses)
#dataFile.close()
#staircase.saveAsPickle(fileName) 

#win.close()
#core.quit()
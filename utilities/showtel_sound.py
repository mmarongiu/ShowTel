from pydub import AudioSegment
from pydub.playback import play
from gtts import gTTS

def voice_create(frase, language):
    voice = gTTS(text=frase, lang=language)
    voice.save("voice.mp3")
    v_out = AudioSegment.from_mp3("voice.mp3")

    return v_out


def sound_sys():
    alarm0 = AudioSegment.from_wav("./sound/warning_alarm.wav")
    v_alarm = alarm0[:2000]                     # first 2 seconds

    voice0_close = gTTS(text="Sun close!", lang='en')
    voice0_close.save("sunclose.mp3")
    v_close = AudioSegment.from_mp3("sunclose.mp3")
    
    voice0_failure = gTTS(text="Failure status!", lang='en')
    voice0_failure.save("failure.mp3")
    v_failure = AudioSegment.from_mp3("failure.mp3")

    voice0_disconnected = gTTS(text="Disconnected!", lang='en')
    voice0_disconnected.save("disconnected.mp3")
    v_disconnected = AudioSegment.from_mp3("disconnected.mp3")

    voice0_connected = gTTS(text="Connected!", lang='en')
    voice0_connected.save("connected.mp3")
    v_connected = AudioSegment.from_mp3("connected.mp3")

    return v_alarm, v_close, v_failure, v_disconnected, v_connected


def sound_alarm(filename, time):
    # time = int in units of seconds
    delta = int(int(time)*1000)
    alarm0 = AudioSegment.from_wav(filename)
    v_alarm = alarm0[:delta]

    return v_alarm


def sound_voice(filename):
    v_out = AudioSegment.from_mp3(filename)

    return v_out

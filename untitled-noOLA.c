#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include "portaudio.h"
#include "sndfile.h"
#include <string.h>
#include <ncurses.h>
#include "inputlib.h"
#include <math.h>
#include "pffft.h"

#ifndef M_PI
#define M_PI	3.14159265358979323846
#endif

#define FRAMES_PER_BUFFER 512
#define MONO 1
#define STEREO 2
#define THRESHOLD_INCREMENT 1
#define ATTENUATION_INCREMENT 0.5
#define INITIAL_THRESHOLD 8
#define INITIAL_ATTENUATION 3
#define INITIAL_ATTACK 50
#define INITIAL_RELEASE 150
#define TC_INCREMENT 5
#define INITIAL_RATIO 4
#define RATIO_INCREMENT 0.5
#define INITIAL_STATE 0

//data struct
typedef struct
{
	float threshold;
	float attenuation;
	float sampleRate;
	float ratio;
	bool toggle;
	SNDFILE *infile1;
	SNDFILE *infile2;
	SF_INFO sfinfo1;
	SF_INFO sfinfo2;
	int attack_frames;
	int release_frames;
	float attack;
	float release;
	float fs_func;
} paData;

//global readcount
int readcount = 0;

//window functions
float win[FRAMES_PER_BUFFER];

//save buffs
float saveprimaryL[FRAMES_PER_BUFFER/2];
float saveprimaryR[FRAMES_PER_BUFFER/2];
float savesecondaryL[FRAMES_PER_BUFFER/2];
float savesecondaryR[FRAMES_PER_BUFFER/2];
float primaryprev[FRAMES_PER_BUFFER];
float secondaryprev[FRAMES_PER_BUFFER];
float magprimary;
float phprimary;
float magsecondary;
float phsecondary;
float magresultant;
float phresultant;
float thetaprimary;
float thetasecondary;
float ratioprimary;
float ratiosecondary;
float primaryproj;
float secondaryproj;
float delta;
float deltaratio;
float newratio;
float newthresh;
float secondaryatten;
int currL[FRAMES_PER_BUFFER/2];
int currR[FRAMES_PER_BUFFER/2];
int prevL[FRAMES_PER_BUFFER/2];
int prevR[FRAMES_PER_BUFFER/2];
int frame_noL[FRAMES_PER_BUFFER/2];
int frame_noR[FRAMES_PER_BUFFER/2];
int currL_sh[FRAMES_PER_BUFFER/2];
int currR_sh[FRAMES_PER_BUFFER/2];
int prevL_sh[FRAMES_PER_BUFFER/2];
int prevR_sh[FRAMES_PER_BUFFER/2];
int frame_noL_sh[FRAMES_PER_BUFFER/2];
int frame_noR_sh[FRAMES_PER_BUFFER/2];
float atkptrL[FRAMES_PER_BUFFER/2][512];
float rlsptrL[FRAMES_PER_BUFFER/2][512];
float atkptrR[FRAMES_PER_BUFFER/2][512];
float rlsptrR[FRAMES_PER_BUFFER/2][512];
int HOP_SIZE;

//callback function declaration
static int paCallback( const void *inputBuffer, void *outputBuffer,
				unsigned long framesPerBuffer,
				const PaStreamCallbackTimeInfo* timeInfo,
				PaStreamCallbackFlags statusFlags,
				void *userData );

void logspace(const float x1, const float x2, const int size, float *array);
void linspace(const float x1, const float x2, const int size, float *array);

void linattack(const float input, const float ratio, const float threshold, const int attack_frames, float *array);
void linrelease(const float input, const float ratio, const float threshold, const int release_frames, float *array);
void logattack(const float input, const float ratio, const float threshold, const int attack_frames, float *array);
void logrelease(const float input, const float ratio, const float threshold, const int release_frames, float *array);
void expattack(const float input, const float ratio, const float threshold, const int attack_frames, float *array);
void exprelease(const float input, const float ratio, const float threshold, const int release_frames, float *array);

void (*attack)(const float, const float, const float, const int, float *);
void (*release)(const float, const float, const float, const int, float *);

int main(int argc, char **argv)
{
	paData data;
	PaStream *stream;
	PaError err;
	PaStreamParameters outputParams;

	//main window
	for (int i = 0; i < FRAMES_PER_BUFFER; i++)
	{
		win[i] = sin((M_PI/(float)FRAMES_PER_BUFFER)*(float)i);
	}

	if (argc != 5)
	{
		printf ("\nPlease input two audio files as input arguments: \n");
		printf ("	Usage : AutoEQ <filename1> <filename2> <attack> <release>\n");
		exit(1);
	}

	//function pointers
	if (strcmp(argv[3], "lin") == 0)
	{
		attack = &linattack;
	}
	else if (strcmp(argv[3], "log") == 0)
	{
		attack = &logattack;
	}
	else if (strcmp(argv[3], "exp") == 0)
	{
		attack = &expattack;
	}
	else{
		printf("Please enter lin, log, or exp for time constant curves.");
		exit(1);
	}

	if (strcmp(argv[4], "lin") == 0)
	{
		release = &linrelease;
	}
	else if (strcmp(argv[4], "log") == 0)
	{
		release = &logrelease;
	}
	else if (strcmp(argv[4], "exp") == 0)
	{
		release = &exprelease;
	}
	else{
		printf("Please enter lin, log, or exp for time constant curves.");
		exit(1);
	}

	//set everything to zero
	memset(&data.sfinfo1, 0, sizeof(data.sfinfo1));
	memset(&data.sfinfo2, 0, sizeof(data.sfinfo2));
	memset(currL,0,sizeof(currL));
	memset(prevL,0,sizeof(prevL));
	memset(currR,0,sizeof(currR));
	memset(prevL,0,sizeof(prevL));
	memset(frame_noL,0,sizeof(frame_noL));
	memset(frame_noR,0,sizeof(frame_noR));
	memset(currL_sh,0,sizeof(currL_sh));
	memset(prevL_sh,0,sizeof(prevL_sh));
	memset(currR_sh,0,sizeof(currR_sh));
	memset(prevL_sh,0,sizeof(prevL_sh));
	memset(frame_noL_sh,0,sizeof(frame_noL_sh));
	memset(frame_noR_sh,0,sizeof(frame_noR_sh));
	//memset(atk,0,sizeof(atk));
	//memset(rls,0,sizeof(rls));

	//open files
	if ((data.infile1 = sf_open(argv[1], SFM_READ, &data.sfinfo1)) == NULL)
	{
		printf("Error: could not open file %s\n", argv[1]);
		puts(sf_strerror(NULL));
		exit(1);
	}

	if ((data.infile2 = sf_open(argv[2], SFM_READ, &data.sfinfo2)) == NULL)
	{
		printf("Error: could not open file %s\n", argv[2]);
		puts(sf_strerror(NULL));
		exit(1);
	}

	printf("Audio file 1: Frames: %d Channels: %d Sample Rate %d\n",
		(int)data.sfinfo1.frames, data.sfinfo1.channels, data.sfinfo1.samplerate);

	printf("Audio file 2: Frames: %d Channels: %d Sample Rate %d\n",
		(int)data.sfinfo2.frames, data.sfinfo2.channels, data.sfinfo2.samplerate);

	//check for length of files
	if (data.sfinfo1.frames != data.sfinfo2.frames)
	{
		printf("Error: audio files must have the same length.\nExiting.\n");
		exit(1);
	}

	//check for channels
	if (data.sfinfo1.channels != data.sfinfo2.channels)
	{
		printf("Error: audio files must have the same number of channels.\nExiting.\n");
		exit(1);
	}

	//check for sample rate
	if (data.sfinfo1.samplerate != data.sfinfo2.samplerate)
	{
		printf("Error: audio files must have the same sample rate.\nExiting.\n");
		exit(1);
	}

	HOP_SIZE = FRAMES_PER_BUFFER/2;
	data.fs_func = (float)HOP_SIZE/(float)data.sfinfo1.samplerate;


	data.sampleRate = data.sfinfo1.samplerate;
	data.threshold = INITIAL_THRESHOLD;
	data.attenuation = INITIAL_ATTENUATION;
	data.attack = INITIAL_ATTACK;
	data.release = INITIAL_RELEASE;
	data.ratio = INITIAL_RATIO;
	data.attack_frames = (int)roundf((data.attack/1000.f)/data.fs_func);
	data.release_frames = (int)roundf((data.release/1000.f)/data.fs_func);
	data.toggle = INITIAL_STATE;

	//initialize and start portaudio
	err = Pa_Initialize();
	if (err != paNoError)
	{
		printf("PortAudio error: %s\n", Pa_GetErrorText(err));
		printf("\nExiting.\n");
		exit(1);
	}

	//output stream parameters
	outputParams.device = Pa_GetDefaultOutputDevice();
	outputParams.channelCount = STEREO;
	outputParams.sampleFormat = paFloat32;
	outputParams.suggestedLatency =
		Pa_GetDeviceInfo(outputParams.device)->defaultLowOutputLatency;
	outputParams.hostApiSpecificStreamInfo = NULL;

	err = Pa_OpenStream(&stream,
			NULL,
			&outputParams,
			data.sampleRate,
			FRAMES_PER_BUFFER/2,
			paNoFlag,
			paCallback,
			&data);

	if (err != paNoError)
	{
		printf("PortAudio error: open stream: %s\n", Pa_GetErrorText(err));
		exit(2);
	}

	//start audio stream
	err = Pa_StartStream(stream);
	if (err != paNoError)
	{
		printf("PortAudio error: start stream: %s\n", Pa_GetErrorText(err));
		exit(3);
	}

	//interactive character input
	initscr();
	cbreak();
	noecho();

	char ch;
	ch = '\0';

	mvprintw(0, 0, "Sensitivity: %.2f Attenuation: %.2f Attack time: %.2f Release time: %.2f Ratio: %.2f On: %d \n \n[k/m] increases/decreases threshold\n" \
		"[a/z] increases/decreases attenuation (dB) \n"\
		"[s/x] increases/decreases release time (ms) \n"\
		"[j/n] increases/decreases attack time (ms) \n"\
		"[d/c] increases/decreases ratio (dB) \n"\
		"[v] to toggle the effect on and off \n"\
		"[q] to quit\n", data.threshold, data.attenuation, data.attack, data.release, data.ratio, data.toggle);

	while (ch != 'q')
	{
		ch = getch();
		switch(ch){
		case 'k':
		{
			data.threshold += THRESHOLD_INCREMENT;
			break;
		}
		case 'm':
		{
			data.threshold -= THRESHOLD_INCREMENT;
			if (data.threshold < 0){
				data.threshold = 0;
			}
			break;
		}
		case 'a':
		{
			data.attenuation += ATTENUATION_INCREMENT;
			break;
		}
		case 'z':
		{
			data.attenuation -= ATTENUATION_INCREMENT;
			if (data.attenuation < 0){
				data.attenuation = 0;
			}
			break;
		}
		case 's':
		{
			data.release += TC_INCREMENT;
			if ((int)roundf((data.release/1000.f)/data.fs_func) < 2)
			{
				data.release_frames = 2;
				break;
			}
			else if ((int)roundf((data.release/1000.f)/data.fs_func) > 512)
			{
				data.release_frames = 512;
				break;
			}
			else
			{
				data.release_frames = (int)roundf((data.release/1000.f)/data.fs_func);
				break;
			}
		}
		case 'x':
		{
			data.release -= TC_INCREMENT;
			if ((int)roundf((data.release/1000.f)/data.fs_func) < 2)
			{
				data.release_frames = 2;
				break;
			}
			else if ((int)roundf((data.release/1000.f)/data.fs_func) > 512)
			{
				data.release_frames = 512;
				break;
			}
			else
			{
				data.release_frames = (int)roundf((data.release/1000.f)/data.fs_func);
				break;
			}
		}
		case 'j':
		{
			data.attack += TC_INCREMENT;
			if ((int)roundf((data.attack/1000.f)/data.fs_func) < 2)
			{
				data.attack_frames = 2;
				break;
			}
			else if ((int)roundf((data.attack/1000.f)/data.fs_func) > 512)
			{
				data.attack_frames = 512;
				break;
			}
			else
			{
				data.attack_frames = (int)roundf((data.attack/1000.f)/data.fs_func);
				break;
			}
		}
		case 'n':
		{
			data.attack -= TC_INCREMENT;
			if ((int)roundf((data.attack/1000.f)/data.fs_func) < 2)
			{
				data.attack_frames = 2;
				break;
			}
			else if ((int)roundf((data.attack/1000.f)/data.fs_func) > 512)
			{
				data.attack_frames = 512;
				break;
			}
			else
			{
				data.attack_frames = (int)roundf((data.attack/1000.f)/data.fs_func);
				break;
			}
		}
		case 'd':
		{
			data.ratio += RATIO_INCREMENT;
			break;
		}
		case 'c':
		{
			data.ratio -= RATIO_INCREMENT;
			if (data.ratio < 1)
			{
				data.ratio = 1;
			}
			break;
		}
		case 'v':
		{
				data.toggle = !data.toggle;
				break;
		}
	}
		mvprintw(0,0,"Sensitivity: %.2f Attenuation: %.2f Attack time: %.2f Release time: %.2f Ratio: %.2f On: %d \n ", data.threshold, data.attenuation, data.attack, data.release, data.ratio, data.toggle);
	}

	//end curses
	endwin();

	//close files
	sf_close(data.infile1);
	sf_close(data.infile2);

	//stop stream
	err = Pa_StopStream(stream);
	if (err != paNoError)
	{
		printf("PortAudio error: stop stream: %s\n", Pa_GetErrorText(err));
	}

	//close stream
	err = Pa_CloseStream(stream);
	if (err != paNoError)
	{
		printf("PortAudio error: close stream: %s\n", Pa_GetErrorText(err));
	}

	//terminate PortAudio
	err = Pa_Terminate();
	if (err != paNoError){
		printf("PortAudio error: terminate: %s\n", Pa_GetErrorText(err));
	}

	return 0;
}

//callback function
static int paCallback( const void *inputBuffer, void *outputBuffer,
				unsigned long framesPerBuffer,
				const PaStreamCallbackTimeInfo* timeInfo,
				PaStreamCallbackFlags statusFlags,
				void *userData)
{
	paData *data = (paData*)userData;
	float *out = (float*)outputBuffer;

	//pointers
	float *primarycurr;
	float *secondarycurr;
	float *realprimaryL;
	float *realprimaryR;
	float *realsecondaryL;
	float *realsecondaryR;
	float *cpxprimaryL;
	float *cpxprimaryR;
	float *cpxsecondaryL;
	float *cpxsecondaryR;
	PFFFT_Setup * setup;
	float amount;
	int index_release;
	int index_attack;

	//attack and release arrays
	//float *atkarrayL = atkptrL[0];
	//float *rlsarrayL = rlsptrL[0];
	//float *atkarrayR = atkprtR[0];
	//float *rlsarrayR = rlsptrR[0];

	//malloc read buffers (stereo)
	primarycurr = (float*)malloc(sizeof(float)*FRAMES_PER_BUFFER);
	secondarycurr = (float*)malloc(sizeof(float)*FRAMES_PER_BUFFER);
	
	//pffft setup
	setup = pffft_new_setup(framesPerBuffer*2,PFFFT_REAL);

	//allocate memory
	int Nbytes = sizeof(float) * framesPerBuffer * 2;

	//alloc real / cpx FFT buffers (mono)
	realprimaryL = (float *)pffft_aligned_malloc(Nbytes);
	realprimaryR = (float *)pffft_aligned_malloc(Nbytes);
	realsecondaryL = (float *)pffft_aligned_malloc(Nbytes);
	realsecondaryR = (float *)pffft_aligned_malloc(Nbytes);
	cpxprimaryL = (float *)pffft_aligned_malloc(Nbytes);
	cpxprimaryR = (float *)pffft_aligned_malloc(Nbytes);
	cpxsecondaryL = (float *)pffft_aligned_malloc(Nbytes);
	cpxsecondaryR = (float *)pffft_aligned_malloc(Nbytes);

	//memset to zero
	memset(primarycurr, 0, sizeof(float)*framesPerBuffer*2);
	memset(secondarycurr, 0, sizeof(float)*framesPerBuffer*2);
	memset(realprimaryL, 0, sizeof(float)*framesPerBuffer*2);
	memset(realprimaryR, 0, sizeof(float)*framesPerBuffer*2);
	memset(realsecondaryL, 0, sizeof(float)*framesPerBuffer*2);
	memset(realsecondaryR, 0, sizeof(float)*framesPerBuffer*2);
	memset(cpxprimaryL, 0, sizeof(float)*framesPerBuffer*2);
	memset(cpxprimaryR, 0, sizeof(float)*framesPerBuffer*2);
	memset(cpxsecondaryL, 0, sizeof(float)*framesPerBuffer*2);
	memset(cpxsecondaryR, 0, sizeof(float)*framesPerBuffer*2);

	//read audio into buffers
	readcount = sf_readf_float(data->infile1, primarycurr, framesPerBuffer);
	readcount = sf_readf_float(data->infile2, secondarycurr, framesPerBuffer);

	//if signal reached the end, rewind
	if (readcount < framesPerBuffer)
	{
		sf_seek(data->infile1, 0, SEEK_SET);
		sf_seek(data->infile2, 0, SEEK_SET);
		sf_readf_float(data->infile1, primarycurr+readcount, framesPerBuffer-readcount);
		sf_readf_float(data->infile2, secondarycurr+readcount, framesPerBuffer-readcount);
	}

	//deinterleave if needed
	if (data->sfinfo1.channels == STEREO)
	{

		//deinterleave and add prev block
		for (int i = 0; i < framesPerBuffer; i++)
		{
	 		realprimaryL[i] = primaryprev[i*2];
	 		realprimaryR[i] = primaryprev[i*2+1];
	 		realsecondaryL[i] = secondaryprev[i*2];
	 		realsecondaryR[i] = secondaryprev[i*2+1];
		}

		//deinterleave and add most recent block
		for (int i = 0; i < framesPerBuffer; i++)
		{
			realprimaryL[i+framesPerBuffer] = primarycurr[i*2];
			realprimaryR[i+framesPerBuffer] = primarycurr[i*2+1];
			realsecondaryL[i+framesPerBuffer] = secondarycurr[i*2];
			realsecondaryR[i+framesPerBuffer] = secondarycurr[i*2+1];
		}

		//copy to prev
		for (int i = 0; i < FRAMES_PER_BUFFER; i++)
		{
			primaryprev[i] = primarycurr[i];
			secondaryprev[i] = secondarycurr[i];
		}

		//window input
		for (int i = 0; i < FRAMES_PER_BUFFER; i++)
		{
	 		realprimaryL[i] = realprimaryL[i] * win[i];
	 		realprimaryR[i] = realprimaryR[i] * win[i];
	 		realsecondaryL[i] = realsecondaryL[i] * win[i];
	 		realsecondaryR[i] = realsecondaryR[i] * win[i];
		}

		//fft transform
		pffft_transform_ordered(setup, realprimaryL, cpxprimaryL, 0, PFFFT_FORWARD);
		pffft_transform_ordered(setup, realprimaryR, cpxprimaryR, 0, PFFFT_FORWARD);
		pffft_transform_ordered(setup, realsecondaryL, cpxsecondaryL, 0, PFFFT_FORWARD);
		pffft_transform_ordered(setup, realsecondaryR, cpxsecondaryR, 0, PFFFT_FORWARD);

		if (data->toggle)
		{

		//left channel autoEQ algo for DC and Nyquist (no imaginary component)
		for (int i = 0; i < 2; i++)
		{
			// if primary bin is above threshold
			if (20.f*log10f(sqrt(powf(cpxprimaryL[2*i],2.f))) > data->threshold)
			{
				//get phases and magnitudes
				magprimary = sqrt(powf(cpxprimaryL[2*i],2.f));
				magsecondary = sqrt(powf(cpxsecondaryL[2*i],2.f));
				phprimary = atan2f(0.f, cpxprimaryL[2*i]);
				phsecondary = atan2f(0.f, cpxsecondaryL[2*i]);

				//wrap to 0 to 2pi
				if (phprimary < 0)
					phprimary += 2.f * M_PI;
				if (phsecondary < 0)
					phsecondary += 2.f * M_PI;

				// if phases are the same or point in opposite directions
				if (fabsf(phprimary - phsecondary) == M_PI || phprimary == phsecondary)
				{
					//if primary doesn't achieve xdB over secondary
					if ((20.f*log10f(magprimary)) - (20.f*log10f(magsecondary)) < data->attenuation)
					{
						//attack and release arrays
						attack(magsecondary,data->ratio,magprimary-powf(10.f,data->attenuation/20.f),data->attack_frames,atkptrL[i]);
						release(magsecondary,data->ratio,magprimary-powf(10.f,data->attenuation/20.f),data->release_frames,rlsptrL[i]);

						//set mode and frame_no if previous frame was inactive
						if (prevL_sh[i] == 0)
						{
							currL_sh[i] = 1;
							frame_noL_sh[i] = 0;
						}
						else if (prevL_sh[i] == 1)
						{
							if (frame_noL_sh[i] == data->attack_frames-1)
							{
								currL_sh[i] = 3;
							}
							else
							{
								frame_noL_sh[i]++;
								currL_sh[i] = 1;
							}
						}
						else if (prevL_sh[i] == 2)
						{
							currL_sh[i] = 1;
							amount = 1;
							index_attack = 0;

							//find nearest point
							for (int p = 0; p < data->attack_frames; ++p)
							{
								if (fabsf(atkptrL[i][p]-rlsptrL[i][frame_noL_sh[i]]) < amount)
								{
									amount = fabsf(atkptrL[i][p]-rlsptrL[i][frame_noL_sh[i]]);
									index_attack = p;
								}
							}

							//find array index
							if ((atkptrL[i][index_attack]>rlsptrL[i][frame_noL_sh[i]]) && (index_attack < data->attack_frames-1))
							{
								frame_noL_sh[i] = index_attack + 1;
							}
							else
							{
								frame_noL_sh[i] = index_attack;
							}
						}
						else if (prevL_sh[i] == 3)
						{
							currL_sh[i] = 3;
							frame_noL_sh[i] = data->attack_frames-1;
						}

						//attenuate secondary
						magsecondary = magsecondary * atkptrL[i][frame_noL_sh[i]];
					
						//back to cartesian
						cpxsecondaryL[2*i] = magsecondary * cosf(phsecondary);
            			cpxsecondaryL[2*i+1] = magsecondary * sinf(phsecondary);
					}

					//if need to check for release
					else if ((20.f*log10f(magprimary) - 20.f*log10f(magsecondary) >= data->attenuation) && (prevL_sh[i] == 1 || prevL_sh[i] == 2 || prevL_sh[i] == 3))
					{
						//attack and release arrays
						//attack(20.f*log10(magsecondary),data->ratio,data->threshold,data->attack_frames,atkarray);
						//release(20.f*log10(magsecondary),data->ratio,data->threshold,data->release_frames,rlsarray);
				
						//if was attack, set to release
						if (prevL_sh[i] == 1)
						{
							currL_sh[i] = 2;
							amount = 1;
							index_release = data->release_frames-1;

							//find nearest point
							for (int n = 0; n < data->release_frames; ++n)
							{
								if (fabsf(atkptrL[i][frame_noL_sh[i]]-rlsptrL[i][n]) < amount)
								{
									amount = fabsf(atkptrL[i][frame_noL_sh[i]]-rlsptrL[i][n]);
									index_release = n;
								}
							}

							//find array index
							if ((atkptrL[i][frame_noL_sh[i]]>rlsptrL[i][index_release]) && (index_release < data->release_frames-1))
							{
								frame_noL_sh[i] = index_release + 1;
							}
							else
							{
								frame_noL_sh[i] = index_release;
							}
						}

						//if was release
						else if (prevL_sh[i] == 2)
						{
							if (frame_noL_sh[i] == data->release_frames - 1)
							{
								currL_sh[i] = 0;
							}
							else
							{
								frame_noL_sh[i]++;
								currL_sh[i] = 2;
							}
						}

						//if was hold
						else if (prevL_sh[i] == 3)
						{
							currL_sh[i] = 2;
							frame_noL_sh[i] = 0;
						}

						//release
						if (currL_sh[i] == 2)
						{
							//attenuate secondary
							magsecondary = magsecondary * rlsptrL[i][frame_noL_sh[i]];
					
							//back to cartesian
							cpxsecondaryL[2*i] = magsecondary * cosf(phsecondary);
            				cpxsecondaryL[2*i+1] = magsecondary * sinf(phsecondary);
						}
					}
				}
			}
		}

	//right channel, DC and Nyquist
	for (int i = 0; i < 2; i++)
	{
		// if primary bin is above threshold
		if (20.f*log10f(sqrt(powf(cpxprimaryR[2*i],2.f))) > data->threshold)
		{
			//get phases and magnitudes
			magprimary = sqrt(powf(cpxprimaryR[2*i],2.f));
			magsecondary = sqrt(powf(cpxsecondaryR[2*i],2.f));
			phprimary = atan2f(0.f, cpxprimaryR[2*i]);
			phsecondary = atan2f(0.f, cpxsecondaryR[2*i]);

			//wrap to 0 to 2pi
			if (phprimary < 0)
				phprimary += 2.f * M_PI;
			if (phsecondary < 0)
				phsecondary += 2.f * M_PI;

			// if phases are the same or point in opposite directions
			if (fabsf(phprimary - phsecondary) == M_PI || phprimary == phsecondary)
			{
				//if primary doesn't achieve xdB over secondary
				if ((20.f*log10f(magprimary)) - (20.f*log10f(magsecondary)) < data->attenuation)
				{
					//attack and release arrays
					attack(magsecondary,data->ratio,magprimary-powf(10.f,data->attenuation/20.f),data->attack_frames,atkptrR[i]);
					release(magsecondary,data->ratio,magprimary-powf(10.f,data->attenuation/20.f),data->release_frames,rlsptrR[i]);

					//set mode and frame_no if previous frame was inactive
					if (prevR_sh[i] == 0)
					{
						currR_sh[i] = 1;
						frame_noR_sh[i] = 0;
					}
					else if (prevR_sh[i] == 1)
					{
						if (frame_noR_sh[i] == data->attack_frames-1)
						{
							currR_sh[i] = 3;
						}
						else
						{
							frame_noR_sh[i]++;
							currR_sh[i] = 1;
						}
					}
					else if (prevR_sh[i] == 2)
					{
						currR_sh[i] = 1;
						amount = 1;
						index_attack = 0;

						//find nearest point
						for (int p = 0; p < data->attack_frames; ++p)
						{
							if (fabsf(atkptrR[i][p]-rlsptrR[i][frame_noR_sh[i]]) < amount)
							{
								amount = fabsf(atkptrR[i][p]-rlsptrR[i][frame_noR_sh[i]]);
								index_attack = p;
							}
						}

						//find array index
						if ((atkptrR[i][index_attack]>rlsptrR[i][frame_noR_sh[i]]) && (index_attack < data->attack_frames-1))
						{
							frame_noR_sh[i] = index_attack + 1;
						}
						else
						{
							frame_noR_sh[i] = index_attack;
						}
					}
					else if (prevR_sh[i] == 3)
					{
						currR_sh[i] = 3;
						frame_noR_sh[i] = data->attack_frames-1;
					}

					//attenuate secondary
					magsecondary = magsecondary * atkptrR[i][frame_noR_sh[i]];
					
					//back to cartesian
					cpxsecondaryR[2*i] = magsecondary * cosf(phsecondary);
            		cpxsecondaryR[2*i+1] = magsecondary * sinf(phsecondary);
				}

			//if need to attenuate
			else if ((20.f*log10f(magprimary) - 20.f*log10f(magsecondary) >= data->attenuation) && (prevR_sh[i] == 1 || prevR_sh[i] == 2 || prevR_sh[i] == 3))
			{
				//attack and release arrays
				//attack(20.f*log10(magsecondary),data->ratio,data->threshold,data->attack_frames,atkarray);
				//release(20.f*log10(magsecondary),data->ratio,data->threshold,data->release_frames,rlsarray);
				
				if (prevR_sh[i] == 1)
				{
					currR_sh[i] = 2;
					amount = 1;
					index_release = data->release_frames-1;

					//find nearest point
					for (int n = 0; n < data->release_frames; ++n)
					{
						if (fabsf(atkptrR[i][frame_noR_sh[i]]-rlsptrR[i][n]) < amount)
						{
							amount = fabsf(atkptrR[i][frame_noR_sh[i]]-rlsptrR[i][n]);
							index_release = n;
						}
					}

					//find array index
					if ((atkptrR[i][frame_noR_sh[i]]>rlsptrR[i][index_release]) && (index_release < data->release_frames-1))
					{
						frame_noR_sh[i] = index_release + 1;
					}
					else
					{
						frame_noR_sh[i] = index_release;
					}
				}
				else if (prevR_sh[i] == 2)
				{
					if (frame_noR_sh[i] == data->release_frames - 1)
					{
						currR_sh[i] = 0;
					}
					else
					{
						frame_noR_sh[i]++;
						currR_sh[i] = 2;
					}
				}
				else if (prevR_sh[i] == 3)
				{
					currR_sh[i] = 2;
					frame_noR_sh[i] = 0;
				}

				//release
				if (currR_sh[i] == 2)
				{
					//attenuate secondary
					magsecondary = magsecondary * rlsptrR[i][frame_noR_sh[i]];
					
					//back to cartesian
					cpxsecondaryR[2*i] = magsecondary * cosf(phsecondary);
            		cpxsecondaryR[2*i+1] = magsecondary * sinf(phsecondary);
				}
			}
		}
	}
}

	//left channel, all other bins
	for (int i = 1; i < framesPerBuffer; i++)
	{
		// if primary bin is above threshold
		if (20.f*log10f(sqrt(powf(cpxprimaryL[2*i],2.f) + powf(cpxprimaryL[2*i+1],2.f))) > data->threshold)
		{
			//get phases and magnitudes
			magprimary = sqrt(powf(cpxprimaryL[2*i],2.f) + powf(cpxprimaryL[2*i+1],2.f));
			magsecondary = sqrt(powf(cpxsecondaryL[2*i],2.f) + powf(cpxsecondaryL[2*i+1],2.f));
			phprimary = atan2f(cpxprimaryL[2*i+1], cpxprimaryL[2*i]);
			phsecondary = atan2f(cpxsecondaryL[2*i+1], cpxsecondaryL[2*i]);

			//wrap to 0 to 2pi
			if (phprimary < 0)
				phprimary += 2.f * M_PI;
			if (phsecondary < 0)
				phsecondary += 2.f * M_PI;

			// if phases are the same or point in opposite directions
			if (fabsf(phprimary - phsecondary) == M_PI || phprimary == phsecondary)
			{
				//if primary doesn't achieve xdB over secondary
				if (20.f*log10f(magprimary) - 20.f*log10f(magsecondary) < data->attenuation)
				{
					//attack and release arrays
					attack(magsecondary,data->ratio,magprimary-powf(10.f,data->attenuation/20.f),data->attack_frames,atkptrL[i]);
					release(magsecondary,data->ratio,magprimary-powf(10.f,data->attenuation/20.f),data->release_frames,rlsptrL[i]);

					//set mode and frame_no if previous frame was inactive
					if (prevL[i] == 0)
					{
						currL[i] = 1;
						frame_noL[i] = 0;
					}
					else if (prevL[i] == 1)
					{
						if (frame_noL[i] == data->attack_frames-1)
						{
							currL[i] = 3;
						}
						else
						{
							frame_noL[i]++;
							currL[i] = 1;
						}
					}
					else if (prevL[i] == 2)
					{
						currL[i] = 1;
						amount = 1;
						index_attack = 0;

						//find nearest point
						for (int p = 0; p < data->attack_frames; ++p)
						{
							if (fabsf(atkptrL[i][p]-rlsptrL[i][frame_noL[i]]) < amount)
							{
								amount = fabsf(atkptrL[i][p]-rlsptrL[i][frame_noL[i]]);
								index_attack = p;
							}
						}

						//find array index
						if ((atkptrL[i][index_attack]>rlsptrL[i][frame_noL[i]]) && (index_attack < data->attack_frames-1))
						{
							frame_noL[i] = index_attack + 1;
						}
						else
						{
							frame_noL[i] = index_attack;
						}
					}
					else if (prevL[i] == 3)
					{
						currL[i] = 3;
						frame_noL[i] = data->attack_frames-1;
					}

					//attenuate secondary
					magsecondary = magsecondary * atkptrL[i][frame_noL[i]];
					
					//back to cartesian
					cpxsecondaryL[2*i] = magsecondary * cosf(phsecondary);
            		cpxsecondaryL[2*i+1] = magsecondary * sinf(phsecondary);
				}

				//if need to release
				else if ((20.f*log10f(magprimary) - 20.f*log10f(magsecondary) >= data->attenuation) && (prevL[i] == 1 || prevL[i] == 2 || prevL[i] == 3))
				{																	
					//attack and release arrays
					//attack(20.f*log10(magsecondary),data->ratio,data->threshold,data->attack_frames,atkarray);
					//release(20.f*log10(magsecondary),data->ratio,data->threshold,data->release_frames,rlsarray);
				
					if (prevL[i] == 1)
					{
						currL[i] = 2;
						amount = 1;
						index_release = data->release_frames-1;

						//find nearest point
						for (int n = 0; n < data->release_frames; ++n)
						{
							if (fabsf(atkptrL[i][frame_noL[i]]-rlsptrL[i][n]) < amount)
							{
								amount = fabsf(atkptrL[i][frame_noL[i]]-rlsptrL[i][n]);
								index_release = n;
							}
						}

						//find array index
						if ((atkptrL[i][frame_noL[i]]>rlsptrL[i][index_release]) && (index_release < data->release_frames-1))
						{
							frame_noL[i] = index_release + 1;
						}
						else
						{
							frame_noL[i] = index_release;
						}
					}
					else if (prevL[i] == 2)
					{
						if (frame_noL[i] == data->release_frames - 1)
						{
							currL_sh[i] = 0;
						}
						else
						{
							frame_noL[i]++;
							currL_sh[i] = 2;
						}
					}
					else if (prevL[i] == 3)
					{
						currL_sh[i] = 2;
						frame_noL[i] = 0;
					}

					//release
					if (currL[i] == 2)
					{
						//attenuate secondary
						magsecondary = magsecondary * rlsptrL[i][frame_noL[i]];
					
						//back to cartesian
						cpxsecondaryL[2*i] = magsecondary * cosf(phsecondary);
            			cpxsecondaryL[2*i+1] = magsecondary * sinf(phsecondary);
					}
				}
			}

			//all other conditions of phase
			else
			{
				//get magnitude and phase of resultant
				magresultant = sqrt(powf(cpxprimaryL[2*i]+cpxsecondaryL[2*i],2.f) + powf(cpxprimaryL[2*i+1]+cpxsecondaryL[2*i+1],2.f));
				phresultant = atan2f(cpxprimaryL[2*i+1]+cpxsecondaryL[2*i+1],cpxprimaryL[2*i]+cpxsecondaryL[2*i]);

				//wrap to 0 to 2pi
				if (phresultant < 0)
					phresultant += 2.f * M_PI;
				if (phresultant < 0)
					phresultant += 2.f * M_PI;

				//calculate delta theta
				thetaprimary = fabsf(phprimary - phresultant);
				thetasecondary = fabsf(phsecondary - phresultant);

				//calculate ratio of scalar projection to hypotenuse
				ratioprimary = fabsf(cosf(thetaprimary));
				ratiosecondary = fabsf(cosf(thetasecondary));

				//calculate scalar projection magnitudes
				primaryproj = magprimary*ratioprimary;
				secondaryproj = magsecondary*ratiosecondary;

				//if scalar projection of primary is less xdB above that of secondary
				if (20.f*log10f(primaryproj) - 20.f*log10f(secondaryproj) < data->attenuation)
				{
					//calculate desired secondary contribution
					secondaryatten = primaryproj - powf(10.f,data->attenuation/20.f);

					//new ratio
					newratio = fabsf(cosf(thetasecondary));

					//secondary threshold
					newthresh = 1.f/(newratio/secondaryatten);

					//attack and release arrays
					attack(20.f*log10f(magsecondary),data->ratio,newthresh,data->attack_frames,atkptrL[i]);
					release(20.f*log10f(magsecondary),data->ratio,newthresh,data->release_frames,rlsptrL[i]);

					//set mode and frame_no if previous frame was inactive
					if (prevL[i] == 0)
					{
						currL[i] = 1;
						frame_noL[i] = 0;
					}
					else if (prevL[i] == 1)
					{
						if (frame_noL[i] == data->attack_frames-1)
						{
							currL[i] = 3;
						}
						else
						{
							frame_noL[i]++;
							currL[i] = 1;
						}
					}
					else if (prevL[i] == 2)
					{
						currL[i] = 1;
						amount = 1;
						index_attack = 0;

						//find nearest point
						for (int p = 0; p < data->attack_frames; ++p)
						{
							if (fabsf(atkptrL[i][p]-rlsptrL[i][frame_noL[i]]) < amount)
							{
								amount = fabsf(atkptrL[i][p]-rlsptrL[i][frame_noL[i]]);
								index_attack = p;
							}
						}

						//find array index
						if ((atkptrL[i][index_attack]>rlsptrL[i][frame_noL[i]]) && (index_attack < data->attack_frames-1))
						{
							frame_noL[i] = index_attack + 1;
						}
						else
						{
							frame_noL[i] = index_attack;
						}
					}
					else if (prevL[i] == 3)
					{
						currL[i] = 3;
						frame_noL[i] = data->attack_frames-1;
					}

					//attenuate secondary
					magsecondary = magsecondary * atkptrL[i][frame_noL[i]];
					
					//back to cartesian
					cpxsecondaryL[2*i] = magsecondary * cosf(phsecondary);
            		cpxsecondaryL[2*i+1] = magsecondary * sinf(phsecondary);
				}

				//if need to release
				else if ((20.f*log10f(primaryproj) - 20.f*log10f(secondaryproj) >= data->attenuation) && (prevL[i] == 1 || prevL[i] == 2 || prevL[i] == 3))
				{
				//attack and release arrays
				//attack(20.f*log10(magsecondary),data->ratio,data->threshold,data->attack_frames,atkarray);
				//release(20.f*log10(magsecondary),data->ratio,data->threshold,data->release_frames,rlsarray);

				if (prevL[i] == 1)
				{
					currL[i] = 2;
					amount = 1;
					index_release = data->release_frames-1;

					//find nearest point
					for (int n = 0; n < data->release_frames; ++n)
					{
						if (fabsf(atkptrL[i][frame_noL[i]]-rlsptrL[i][n]) < amount)
						{
							amount = fabsf(atkptrL[i][frame_noL[i]]-rlsptrL[i][n]);
							index_release = n;
						}
					}

					//find array index
					if ((atkptrL[i][frame_noL[i]]>rlsptrL[i][index_release]) && (index_release < data->release_frames-1))
					{
						frame_noL[i] = index_release + 1;
					}
					else
					{
						frame_noL[i] = index_release;
					}
				}
				else if (prevL[i] == 2)
				{
					if (frame_noL[i] == data->release_frames - 1)
					{
						currL[i] = 0;
					}
					else
					{
						frame_noL[i]++;
						currL[i] = 2;
					}
				}
				else if (prevL[i] == 3)
				{
					currL[i] = 2;
					frame_noL[i] = 0;
				}

				//release
				if (currL[i] == 2)
				{
					//attenuate secondary
					magsecondary = magsecondary * rlsptrL[i][frame_noL[i]];
					
					//back to cartesian
					cpxsecondaryL[2*i] = magsecondary * cosf(phsecondary);
            		cpxsecondaryL[2*i+1] = magsecondary * sinf(phsecondary);
				}
			}
		}
	}
}

	//right channel, all other bins
	for (int i = 1; i < framesPerBuffer; i++)
	{
		// if primary bin is above threshold
		if (20.f*log10f(sqrt(powf(cpxprimaryR[2*i],2.f) + powf(cpxprimaryR[2*i+1],2.f))) > data->threshold)
		{
			//get phases and magnitudes
			magprimary = sqrt(powf(cpxprimaryR[2*i],2.f) + powf(cpxprimaryR[2*i+1],2.f));
			magsecondary = sqrt(powf(cpxsecondaryR[2*i],2.f) + powf(cpxsecondaryR[2*i+1],2.f));
			phprimary = atan2f(cpxprimaryR[2*i+1], cpxprimaryR[2*i]);
			phsecondary = atan2f(cpxsecondaryR[2*i+1], cpxsecondaryR[2*i]);

			//wrap to 0 to 2pi
			if (phprimary < 0)
				phprimary += 2.f * M_PI;
			if (phsecondary < 0)
				phsecondary += 2.f * M_PI;

			// if phases are the same or point in opposite directions
			if (fabsf(phprimary - phsecondary) == M_PI || phprimary == phsecondary)
			{
				//if primary doesn't achieve xdB over secondary
				if ((20.f*log10f(magprimary)) - (20.f*log10f(magsecondary)) < data->attenuation)
				{
					//attack and release arrays
					attack(20.f*log10f(magsecondary),data->ratio,magprimary-powf(10.f,data->attenuation/20.f),data->attack_frames,atkptrR[i]);
					release(20.f*log10f(magsecondary),data->ratio,magprimary-powf(10.f,data->attenuation/20.f),data->release_frames,rlsptrR[i]);

					//set mode and frame_no if previous frame was inactive
					if (prevR[i] == 0)
					{
						currR[i] = 1;
						frame_noR[i] = 0;
					}
					else if (prevR[i] == 1)
					{
						if (frame_noR[i] == data->attack_frames-1)
						{
							currR[i] = 3;
						}
						else
						{
							frame_noR[i]++;
							currR[i] = 1;
						}
					}
					else if (prevR[i] == 2)
					{
						currR[i] = 1;
						amount = 1;
						index_attack = 0;

						//find nearest point
						for (int p = 0; p < data->attack_frames; ++p)
						{
							if (fabsf(atkptrR[i][p]-rlsptrR[i][frame_noR[i]]) < amount)
							{
								amount = fabsf(atkptrR[i][p]-rlsptrR[i][frame_noR[i]]);
								index_attack = p;
							}
						}

						//find array index
						if ((atkptrR[i][index_attack]>rlsptrR[i][frame_noR[i]]) && (index_attack < data->attack_frames-1))
						{
							frame_noR[i] = index_attack + 1;
						}
						else
						{
							frame_noR[i] = index_attack;
						}
					}
					else if (prevR[i] == 3)
					{
						currR[i] = 3;
						frame_noR[i] = data->attack_frames-1;
					}

					//attenuate secondary
					magsecondary = magsecondary * atkptrR[i][frame_noR[i]];
					
					//back to cartesian
					cpxsecondaryR[2*i] = magsecondary * cosf(phsecondary);
            		cpxsecondaryR[2*i+1] = magsecondary * sinf(phsecondary);
				}
			}
			else if ((20.f*log10f(magprimary) - 20.f*log10f(magsecondary) >= data->attenuation) && (prevR[i] == 1 || prevR[i] == 2 || prevR[i] == 3))
			{
				//attack and release arrays
				//attack(20.f*log10(magsecondary),data->ratio,data->threshold,data->attack_frames,atkarray);
				//release(20.f*log10(magsecondary),data->ratio,data->threshold,data->release_frames,rlsarray);
				
				if (prevR[i] == 1)
				{
					currR[i] = 2;
					amount = 1;
					index_release = data->release_frames-1;

					//find nearest point
					for (int n = 0; n < data->release_frames; ++n)
					{
						if (fabsf(atkptrR[i][frame_noR[i]]-rlsptrR[i][n]) < amount)
						{
							amount = fabsf(atkptrR[i][frame_noR[i]]-rlsptrR[i][n]);
							index_release = n;
						}
					}

					//find array index
					if ((atkptrR[i][frame_noR[i]]>rlsptrR[i][index_release]) && (index_release < data->release_frames-1))
					{
						frame_noR[i] = index_release + 1;
					}
					else
					{
						frame_noR[i] = index_release;
					}
				}
				else if (prevR[i] == 2)
				{
					if (frame_noR[i] == data->release_frames - 1)
					{
						currR[i] = 0;
					}
					else
					{
						frame_noR[i]++;
						currR[i] = 2;
					}
				}
				else if (prevR[i] == 3)
				{
					currR[i] = 2;
					frame_noR[i] = 0;
				}

				//release
				if (currR[i] == 2)
				{
					//attenuate secondary
					magsecondary = magsecondary * rlsptrR[i][frame_noR[i]];
					
					//back to cartesian
					cpxsecondaryR[2*i] = magsecondary * cosf(phsecondary);
            		cpxsecondaryR[2*i+1] = magsecondary * sinf(phsecondary);
				}
			}

			//all other conditions of phase
			else
			{
				//get magnitude and phase of resultant
				magresultant = sqrt(powf(cpxprimaryR[2*i]+cpxsecondaryR[2*i],2.f) + powf(cpxprimaryR[2*i+1]+cpxsecondaryR[2*i+1],2.f));
				phresultant = atan2f(cpxprimaryR[2*i+1]+cpxsecondaryR[2*i+1],cpxprimaryR[2*i]+cpxsecondaryR[2*i]);

				//wrap to 0 to 2pi
				if (phresultant < 0)
					phresultant += 2.f * M_PI;
				if (phresultant < 0)
					phresultant += 2.f * M_PI;

				//calculate delta theta
				thetaprimary = fabsf(phprimary - phresultant);
				thetasecondary = fabsf(phsecondary - phresultant);

				//calculate ratio of scalar projection to hypotenuse
				ratioprimary = fabsf(cosf(thetaprimary));
				ratiosecondary = fabsf(cosf(thetasecondary));

				//calculate scalar projection magnitudes
				primaryproj = magprimary*ratioprimary;
				secondaryproj = magsecondary*ratiosecondary;

				//if scalar projection of primary is less xdB above that of secondary
				if (20.f*log10f(primaryproj) - 20.f*log10f(secondaryproj) < data->attenuation)
				{

					//calculate desired secondary contribution
					secondaryatten = primaryproj - powf(10.f,data->attenuation/20.f);

					//new ratio
					newratio = fabsf(cosf(thetasecondary));

					//secondary threshold
					newthresh = 1.f/(newratio/secondaryatten);

					//attack and release arrays
					attack(20.f*log10f(magsecondary),data->ratio,newthresh,data->attack_frames,atkptrR[i]);
					release(20.f*log10f(magsecondary),data->ratio,newthresh,data->release_frames,rlsptrR[i]);

					//set mode and frame_no if previous frame was inactive
					if (prevR[i] == 0)
					{
						currR[i] = 1;
						frame_noR[i] = 0;
					}
					else if (prevR[i] == 1)
					{
						if (frame_noR[i] == data->attack_frames-1)
						{
							currR[i] = 3;
						}
						else
						{
							frame_noL[i]++;
							currR[i] = 1;
						}
					}
					else if (prevR[i] == 2)
					{
						currR[i] = 1;
						amount = 1;
						index_attack = 0;

						//find nearest point
						for (int p = 0; p < data->attack_frames; ++p)
						{
							if (fabsf(atkptrR[i][p]-rlsptrR[i][frame_noL[i]]) < amount)
							{
								amount = fabsf(atkptrR[i][p]-rlsptrR[i][frame_noR[i]]);
								index_attack = p;
							}
						}

						//find array index
						if ((atkptrR[i][index_attack]>rlsptrR[i][frame_noL[i]]) && (index_attack < data->attack_frames-1))
						{
							frame_noR[i] = index_attack + 1;
						}
						else
						{
							frame_noR[i] = index_attack;
						}
					}
					else if (prevR[i] == 3)
					{
						currR[i] = 3;
						frame_noR[i] = data->attack_frames-1;
					}

					//attenuate secondary
					magsecondary = magsecondary * atkptrR[i][frame_noR[i]];
					
					//back to cartesian
					cpxsecondaryR[2*i] = magsecondary * cosf(phsecondary);
            		cpxsecondaryR[2*i+1] = magsecondary * sinf(phsecondary);
				}
			else if ((20.f*log10f(primaryproj) - 20.f*log10f(secondaryproj) >= data->attenuation) && (prevR[i] == 1 || prevR[i] == 2 || prevR[i] == 3))
			{
				//attack and release arrays
				//attack(20.f*log10(magsecondary),data->ratio,data->threshold,data->attack_frames,atkarray);
				//release(20.f*log10(magsecondary),data->ratio,data->threshold,data->release_frames,rlsarray);
				
				if (prevR[i] == 1)
				{
					currR[i] = 2;
					amount = 1;
					index_release = data->release_frames-1;

					//find nearest point
					for (int n = 0; n < data->release_frames; ++n)
					{
						if (fabsf(atkptrR[i][frame_noR[i]]-rlsptrR[i][n]) < amount)
						{
							amount = fabsf(atkptrR[i][frame_noR[i]]-rlsptrR[i][n]);
							index_release = n;
						}
					}

					//find array index
					if ((atkptrR[i][frame_noR[i]]>rlsptrR[i][index_release]) && (index_release < data->release_frames-1))
					{
						frame_noR[i] = index_release + 1;
					}
					else
					{
						frame_noR[i] = index_release;
					}
				}
				else if (prevR[i] == 2)
				{
					if (frame_noL[i] == data->release_frames - 1)
					{
						currR[i] = 0;
					}
					else
					{
						frame_noR[i]++;
						currR[i] = 2;
					}
				}
				else if (prevR[i] == 3)
				{
					currR[i] = 2;
					frame_noR[i] = 0;
				}

				//release
				if (currR[i] == 2)
				{
					//attenuate secondary
					magsecondary = magsecondary * rlsptrR[i][frame_noR[i]];
					
					//back to cartesian
					cpxsecondaryR[2*i] = magsecondary * cosf(phsecondary);
            		cpxsecondaryR[2*i+1] = magsecondary * sinf(phsecondary);
				}
			}
		}
	}
}
}
	//ifft transform
	pffft_transform_ordered(setup, cpxprimaryL, realprimaryL, 0, PFFFT_BACKWARD);
	pffft_transform_ordered(setup, cpxprimaryR, realprimaryR, 0, PFFFT_BACKWARD);
	pffft_transform_ordered(setup, cpxsecondaryL, realsecondaryL, 0, PFFFT_BACKWARD);
	pffft_transform_ordered(setup, cpxsecondaryR, realsecondaryR, 0, PFFFT_BACKWARD);

	//window output
	for (int i = 0; i < FRAMES_PER_BUFFER; i++)
	{
	 	realprimaryL[i] = realprimaryL[i] * win[i];
	 	realprimaryR[i] = realprimaryR[i] * win[i];
	 	realsecondaryL[i] = realsecondaryL[i] * win[i];
	 	realsecondaryR[i] = realsecondaryR[i] * win[i];
	}

	//normalize and copy to output
	for (int i = 0; i < framesPerBuffer; i++)
	{
		out[2*i] = (realprimaryL[i] + realsecondaryL[i] + saveprimaryL[i] + savesecondaryL[i]) / (float)FRAMES_PER_BUFFER;
		out[2*i+1] = (realprimaryR[i] + realsecondaryR[i] + saveprimaryR[i] + savesecondaryR[i]) / (float)FRAMES_PER_BUFFER;
	}

	//save 2nd half
	for (int i = framesPerBuffer; i < FRAMES_PER_BUFFER; i++)
	{
		saveprimaryL[i-framesPerBuffer] = realprimaryL[i];
		saveprimaryR[i-framesPerBuffer] = realprimaryR[i];
		savesecondaryL[i-framesPerBuffer] = realsecondaryL[i];
		savesecondaryR[i-framesPerBuffer] = realsecondaryR[i];
	}

	//update modes
	for (int i = 0; i < framesPerBuffer; ++i)
	{
		prevL[i] = currL[i];
		prevR[i] = currR[i];
		prevL_sh[i] = currL_sh[i];
		prevR_sh[i] = currR_sh[i];
	}

	//clean up
	pffft_destroy_setup(setup);
	pffft_aligned_free(realprimaryL);
	pffft_aligned_free(realprimaryR);
	pffft_aligned_free(realsecondaryL);
	pffft_aligned_free(realsecondaryR);
	pffft_aligned_free(cpxprimaryL);
	pffft_aligned_free(cpxprimaryR);
	pffft_aligned_free(cpxsecondaryL);
	pffft_aligned_free(cpxsecondaryR);
	free(primarycurr);
	free(secondarycurr);

	}
	return 0;
}

void logspace(const float x1, const float x2, const int size, float *array)
{
		if(size == 1)
		{
			array[0] = x2;
		}
		else{
		float delta_x = (x2-x1)/(size-1);

		for(int ii = 0; ii < size; ii++)
		{
			array[ii] = powf(10.0,(x1 + delta_x*ii));
		}
		}
		return;
}

void linspace(const float x1, const float x2, const int size, float *array)
{
		if(size == 1)
		{
			array[0] = x2;
		}
		else{
		float delta_x = (x2-x1)/(size-1);

		for(int ii = 0; ii < size; ii++)
		{
			array[ii] = x1 + ii*delta_x;
		}

		array[size-1] = x2;
		}
		return;
}

void linattack(const float input, const float ratio, const float threshold, const int attack_frames, float *array)
{
	linspace(1.f,(threshold+(input-threshold)/ratio)/input,attack_frames, array);
	return;
}

void linrelease(const float input, const float ratio, const float threshold, const int release_frames, float *array)
{
	linspace((threshold+(input-threshold)/ratio)/input,1.f,release_frames, array);
	return;
}

void logattack(const float input, const float ratio, const float threshold, const int attack_frames, float *array)
{
	logspace(log10f((threshold+(input-threshold)/ratio)/input),0,attack_frames, array);

	float factor = (1+(threshold+(input-threshold)/ratio)/input);

	for (int i = 0; i < attack_frames; ++i)
	{
		array[i] = factor - array[i];
	}
	return;
}

void logrelease(const float input, const float ratio, const float threshold, const int release_frames, float *array)
{
	logspace(log10f((threshold+(input-threshold)/ratio)/input),0,release_frames, array);return;
}

void expattack(const float input, const float ratio, const float threshold, const int attack_frames, float *array)
{
	logspace(0,log10f((threshold+(input-threshold)/ratio)/input),attack_frames, array);return;
}

void exprelease(const float input, const float ratio, const float threshold, const int release_frames, float *array)
{
	logspace(0,log10f((threshold+(input-threshold)/ratio)/input),release_frames, array);

	float factor = (1+(threshold+(input-threshold)/ratio)/input);

	for (int i = 0; i < release_frames; ++i)
	{
		array[i] = factor - array[i];
	}
	return;
}
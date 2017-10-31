//Eric Zhang - Automatic Surgical EQ
//Portaudio code by M. Farbood

#include <stdlib.h>
#include <stdio.h>
#include <portaudio.h>
#include <sndfile.h>
#include <string.h>
#include <ncurses.h>
#include "inputlib.h"
#include <math.h>
#include "fftw3.h"

#define FRAMES_PER_BUFFER 512
#define MONO 1
#define STEREO 2
#define SENSITIVITY_INCREMENT 1
#define ATTENUATION_INCREMENT 0.5
#define INITIAL_SENSITIVITY 5
#define INITIAL_ATTENUATION 3

//data struct
typedef struct
{
	float sensitivity;
	float attenuation;
	float sampleRate;
	SNDFILE *infile1;
	SNDFILE *infile2;
	SF_INFO sfinfo1;
	SF_INFO sfinfo2;
} paData;

//global variables
unsigned flags = FFTW_MEASURE;

//wisdom file
const char *filename = "surgEQ.wis";

//check for the start of the signal to take a half frame to window
bool start = TRUE;

//initialize sample offset counter for hop size rewinding
int framesleft = FRAMES_PER_BUFFER;

//global readcount
int readcount;

//hop size
int hop = FRAMES_PER_BUFFER/2;

//window functions
float win[FRAMES_PER_BUFFER];
float coswin[FRAMES_PER_BUFFER/2];

//callback function declaration
static int paCallback( const void *inputBuffer, void *outputBuffer,
				unsigned long framesPerBuffer,
				const PaStreamCallbackTimeInfo* timeInfo,
				PaStreamCallbackFlags statusFlags,
				void *userData );

//global pointers for prev arrays for OLA
fftw_complex *primaryprevL = NULL;
fftw_complex *primaryprevR = NULL;
fftw_complex *secondaryprevL = NULL;
fftw_complex *secondaryprevR = NULL;

int main( int argc, char **argv)
{
	paData data;
	PaStream *stream;
	PaError err;
	PaStreamParameters outputParams;

	//set windows to zero
	memset(win,0,sizeof(float)*FRAMES_PER_BUFFER);
	memset(coswin,0,sizeof(float)*FRAMES_PER_BUFFER/2);

	//main window
	for (int i = 0; i < FRAMES_PER_BUFFER; i++)
	{
		win[i] = sin((M_PI/(float)FRAMES_PER_BUFFER)*(float)i);
	}

	//half cosine window
	for (int i = FRAMES_PER_BUFFER/2; i < FRAMES_PER_BUFFER; i++)
	{
		coswin[i] = win[i];
	}

	if (argc != 3)
	{
		printf ("\nPlease input two audio files as input arguments: \n");
		printf ("	Usage : AutoEQ <filename1> <filename2>\n");
		exit(1);
	}

	//set everything to zero
	memset(&data.sfinfo1, 0, sizeof(data.sfinfo1));
	memset(&data.sfinfo2, 0, sizeof(data.sfinfo2));

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

	data.sampleRate = data.sfinfo1.samplerate;
	data.sensitivity = INITIAL_SENSITIVITY;
	data.attenuation = 20 * log10(INITIAL_ATTENUATION);

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

	mvprintw(0, 0, "Sensitivity: %d Attenuation: %.2f \n \n[a/z] increases/decreases sensitivity (dB) \n" \
		"[k/m] increases/decreases attenuation (dB) \n"\
		"[q] to quit\n", data.sensitivity, data.attenuation);

	while (ch != 'q')
	{
		ch = getch();
		switch(ch){
		case 'k':
		{
			data.sensitivity += SENSITIVITY_INCREMENT;
			break;
		}
		case 'm':
		{
			data.sensitivity -= SENSITIVITY_INCREMENT;
			if (data.sensitivity < 0){
				data.sensitivity = 0;
			}
			break;
		}
		case 'a':
		{
			data.attenuation += 20 * log10(ATTENUATION_INCREMENT);
			break;
		}
		case 'z':
		{
			data.attenuation -= 20 * log10(ATTENUATION_INCREMENT);
			if (data.attenuation < 0)
			{
				data.attenuation = 0;
			}
			break;
		}
	}
		mvprintw(0, 0, "Sensitivity: %d Attenuation: %.2f \n ", data.sensitivity, data.attenuation);
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

	int i;

	//deinterleave if needed
	if (data->sfinfo1.channels == STEREO)
	{

	//malloc complex arrays for fftw
	fftw_complex *cbuffprimaryL = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * FRAMES_PER_BUFFER);
	fftw_complex *cbuffprimaryR = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * FRAMES_PER_BUFFER);
	fftw_complex *outcbuffprimaryL = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * FRAMES_PER_BUFFER);
	fftw_complex *outcbuffprimaryR = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * FRAMES_PER_BUFFER);

	fftw_complex *cbuffsecondaryL = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * FRAMES_PER_BUFFER);
	fftw_complex *cbuffsecondaryR = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * FRAMES_PER_BUFFER);
	fftw_complex *outcbuffsecondaryL = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * FRAMES_PER_BUFFER);
	fftw_complex *outcbuffsecondaryR = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * FRAMES_PER_BUFFER);

	//malloc float arrays for reading
	float *realprimarycurr = (float *)fftw_malloc(sizeof(float)*FRAMES_PER_BUFFER*2);
	float *realprimarycurrL = (float *)fftw_malloc(sizeof(float)*FRAMES_PER_BUFFER);
	float *realprimarycurrR = (float *)fftw_malloc(sizeof(float)*FRAMES_PER_BUFFER);
	float *realsecondarycurr = (float *)fftw_malloc(sizeof(float)*FRAMES_PER_BUFFER*2);
	float *realsecondarycurrL = (float *)fftw_malloc(sizeof(float)*FRAMES_PER_BUFFER);
	float *realsecondarycurrR = (float *)fftw_malloc(sizeof(float)*FRAMES_PER_BUFFER);

	//starting 1/2 block
	float *instartprimary = (float *)fftw_malloc(sizeof(float)*FRAMES_PER_BUFFER);
	float *instartprimaryL = (float *)fftw_malloc(sizeof(float)*FRAMES_PER_BUFFER/2);
	float *instartprimaryR = (float *)fftw_malloc(sizeof(float)*FRAMES_PER_BUFFER/2);
	float *instartsecondary = (float *)fftw_malloc(sizeof(float)*FRAMES_PER_BUFFER);
	float *instartsecondaryL = (float *)fftw_malloc(sizeof(float)*FRAMES_PER_BUFFER/2);
	float *instartsecondaryR = (float *)fftw_malloc(sizeof(float)*FRAMES_PER_BUFFER/2);

	fftw_complex *cstartprimaryL = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * FRAMES_PER_BUFFER/2);
	fftw_complex *cstartprimaryR = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * FRAMES_PER_BUFFER/2);
	fftw_complex *outcstartprimaryL = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * FRAMES_PER_BUFFER/2);
	fftw_complex *outcstartprimaryR = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * FRAMES_PER_BUFFER/2);
	fftw_complex *cstartsecondaryL = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * FRAMES_PER_BUFFER/2);
	fftw_complex *cstartsecondaryR = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * FRAMES_PER_BUFFER/2);
	fftw_complex *outcstartsecondaryL = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * FRAMES_PER_BUFFER/2);
	fftw_complex *outcstartsecondaryR = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * FRAMES_PER_BUFFER/2);

	//malloc complex arrays for saving
	fftw_complex *savecbuffprimaryL = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * FRAMES_PER_BUFFER);
	fftw_complex *savecbuffprimaryR = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * FRAMES_PER_BUFFER);
	fftw_complex *savecbuffsecondaryL = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * FRAMES_PER_BUFFER);
	fftw_complex *savecbuffsecondaryR = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * FRAMES_PER_BUFFER);

	//declare input plans
	fftw_plan planprimarycurrL;
	fftw_plan planprimarycurrR;
	fftw_plan plansecondarycurrL;
	fftw_plan plansecondarycurrR;

	//starting block plans
	fftw_plan planstartprimaryL;
	fftw_plan planstartprimaryR;
	fftw_plan planstartsecondaryL;
	fftw_plan planstartsecondaryR;
	fftw_plan outplanstartprimaryL;
	fftw_plan outplanstartprimaryR;
	fftw_plan outplanstartsecondaryL;
	fftw_plan outplanstartsecondaryR;

	//declare output plans
	fftw_plan outplanprimarycurrL;
	fftw_plan outplanprimarycurrR;
	fftw_plan outplansecondarycurrL;
	fftw_plan outplansecondarycurrR;

	//reset readcount
 	readcount = 0;

	//read audio into buffers
	if (start == TRUE)
	{
		sf_readf_float(data->infile1, instartprimary, FRAMES_PER_BUFFER/2);
		sf_readf_float(data->infile2, instartsecondary, FRAMES_PER_BUFFER/2);
		sf_seek(data->infile1, 0, SEEK_SET);
		sf_seek(data->infile2, 0, SEEK_SET);
		readcount = sf_readf_float(data->infile1, realprimarycurr, FRAMES_PER_BUFFER);
		readcount = sf_readf_float(data->infile2, realsecondarycurr, FRAMES_PER_BUFFER);
	}
	else if (framesleft < hop && start == FALSE)
	{
		sf_seek(data->infile1, -(hop-framesleft), SEEK_END);
		sf_seek(data->infile2, -(hop-framesleft), SEEK_END);
		sf_readf_float(data->infile1, realprimarycurr, hop-framesleft);
		sf_readf_float(data->infile2, realsecondarycurr, hop-framesleft);
		sf_seek(data->infile1, 0, SEEK_SET);
		sf_seek(data->infile2, 0, SEEK_SET);
		sf_readf_float(data->infile1, realprimarycurr+(hop-framesleft), FRAMES_PER_BUFFER-(hop-framesleft));
		sf_readf_float(data->infile2, realsecondarycurr+(hop-framesleft), FRAMES_PER_BUFFER-(hop-framesleft));
	}
	else
	{
		sf_seek(data->infile1, -hop, SEEK_CUR);
		sf_seek(data->infile2, -hop, SEEK_CUR);
		readcount = sf_readf_float(data->infile1, realprimarycurr, FRAMES_PER_BUFFER);
		readcount = sf_readf_float(data->infile2, realsecondarycurr, FRAMES_PER_BUFFER);
	}

	//if signal reached the end, rewind
	if (readcount < FRAMES_PER_BUFFER)
	{
		sf_seek(data->infile1, 0, SEEK_SET);
		sf_seek(data->infile2, 0, SEEK_SET);
		framesleft = sf_readf_float(data->infile1, realprimarycurr+readcount, FRAMES_PER_BUFFER-readcount);
		framesleft = sf_readf_float(data->infile2, realsecondarycurr+readcount, FRAMES_PER_BUFFER-readcount);
	}

	//deinterleave
	if (start == TRUE)
	{
		for (i = 0; i < FRAMES_PER_BUFFER/2; i++)
		{
			instartprimaryL[i] = instartprimary[i*2];
			instartprimaryR[i] = instartprimary[i*2+1];
			instartsecondaryL[i] = instartsecondary[i*2];
			instartsecondaryR[i] = instartsecondary[i*2+1];
		}
		for (i = 0; i < FRAMES_PER_BUFFER; i++)
		{
			realprimarycurrL[i] = realprimarycurr[i*2];
			realprimarycurrR[i] = realprimarycurr[i*2+1];
			realsecondarycurrL[i] = realsecondarycurr[i*2];
			realsecondarycurrR[i] = realsecondarycurr[i*2+1];
		}
	}
	else
	{
		for (i = 0; i < FRAMES_PER_BUFFER; i++)
		{
			realprimarycurrL[i] = realprimarycurr[i*2];
			realprimarycurrR[i] = realprimarycurr[i*2+1];
			realsecondarycurrL[i] = realsecondarycurr[i*2];
			realsecondarycurrR[i] = realsecondarycurr[i*2+1];
		}
	}

	// for (i = 0; i < FRAMES_PER_BUFFER/2; i++)
	// {
	// 	out[i*2] = realprimarycurrL[i] + realsecondarycurrL[i];
	// 	out[i*2+1] = realsecondarycurrR[i] + realsecondarycurrR[i];
	// }

	//window real input and copy to complex arrays
	if (start == TRUE)
	{
		for (i = 0; i < FRAMES_PER_BUFFER; i++)
		{
			cbuffprimaryL[i][0] = win[i] * realprimarycurrL[i];
			cbuffprimaryL[i][1] = 0.f;

			cbuffprimaryR[i][0] = win[i] * realprimarycurrR[i];
			cbuffprimaryR[i][1] = 0.f;

			cbuffsecondaryL[i][0] = win[i] * realsecondarycurrL[i];
			cbuffsecondaryL[i][1] = 0.f;

			cbuffsecondaryR[i][0] = win[i] * realsecondarycurrR[i];
			cbuffsecondaryR[i][1] = 0.f;
		}

		for (i = 0; i < FRAMES_PER_BUFFER/2; i++)
		{
			cstartprimaryL[i][0] = coswin[i] * instartprimaryL[i];
			cstartprimaryL[i][1] = 0.f;

			cstartprimaryR[i][0] = coswin[i] * instartprimaryR[i];
			cstartprimaryR[i][1] = 0.f;

			cstartsecondaryL[i][0] = coswin[i] * instartsecondaryL[i];
			cstartsecondaryL[i][1] = 0.f;

			cstartsecondaryR[i][0] = coswin[i] * instartsecondaryR[i];
			cstartsecondaryR[i][1] = 0.f;
		}
	}
	else
	{
		for (i = 0; i < FRAMES_PER_BUFFER; i++)
		{
			cbuffprimaryL[i][0] = win[i] * realprimarycurrL[i];
			cbuffprimaryL[i][1] = 0.f;

			cbuffprimaryR[i][0] = win[i] * realprimarycurrR[i];
			cbuffprimaryR[i][1] = 0.f;

			cbuffsecondaryL[i][0] = win[i] * realsecondarycurrL[i];
			cbuffsecondaryL[i][1] = 0.f;

			cbuffsecondaryR[i][0] = win[i] * realsecondarycurrR[i];
			cbuffsecondaryR[i][1] = 0.f;
		}
	}

	// //import wisdom file
	if (fftw_import_wisdom_from_filename(filename) == 0)
	{
		printf("Error importing wisdom file.");
		exit(1);
	}

	//create forward plans
	if (start == TRUE)
	{
		planprimarycurrL = fftw_plan_dft_1d(FRAMES_PER_BUFFER, cbuffprimaryL, outcbuffprimaryL, FFTW_FORWARD, flags);
		planprimarycurrR = fftw_plan_dft_1d(FRAMES_PER_BUFFER, cbuffprimaryR, outcbuffprimaryR, FFTW_FORWARD, flags);
		plansecondarycurrL = fftw_plan_dft_1d(FRAMES_PER_BUFFER, cbuffsecondaryL, outcbuffsecondaryL, FFTW_FORWARD, flags);
		plansecondarycurrR = fftw_plan_dft_1d(FRAMES_PER_BUFFER, cbuffsecondaryR, outcbuffsecondaryR, FFTW_FORWARD, flags);
		planstartprimaryL = fftw_plan_dft_1d(FRAMES_PER_BUFFER/2, cstartprimaryL, outcstartprimaryL, FFTW_FORWARD, flags);
		planstartprimaryR = fftw_plan_dft_1d(FRAMES_PER_BUFFER/2, cstartprimaryR, outcstartprimaryR, FFTW_FORWARD, flags);
		planstartsecondaryL = fftw_plan_dft_1d(FRAMES_PER_BUFFER/2, cstartsecondaryL, outcstartsecondaryL, FFTW_FORWARD, flags);
		planstartsecondaryR = fftw_plan_dft_1d(FRAMES_PER_BUFFER/2, cstartsecondaryR, outcstartsecondaryR, FFTW_FORWARD, flags);
	}
	else
	{
		planprimarycurrL = fftw_plan_dft_1d(FRAMES_PER_BUFFER, cbuffprimaryL, outcbuffprimaryL, FFTW_FORWARD, flags);
		planprimarycurrR = fftw_plan_dft_1d(FRAMES_PER_BUFFER, cbuffprimaryR, outcbuffprimaryR, FFTW_FORWARD, flags);
		plansecondarycurrL = fftw_plan_dft_1d(FRAMES_PER_BUFFER, cbuffsecondaryL, outcbuffsecondaryL, FFTW_FORWARD, flags);
		plansecondarycurrR = fftw_plan_dft_1d(FRAMES_PER_BUFFER, cbuffsecondaryR, outcbuffsecondaryR, FFTW_FORWARD, flags);
	}

	//export wisdom file
	// if (fftw_export_wisdom_to_filename(filename) == 0)
	// {
	// 	printf("Error exporting wisdom file.");
	// 	exit(1);
	// }

	//execute forward plans
	if (start == TRUE)
	{
		fftw_execute(planprimarycurrL);
		fftw_execute(planprimarycurrR);
		fftw_execute(plansecondarycurrL);
		fftw_execute(plansecondarycurrR);
		fftw_execute(planstartprimaryL);
		fftw_execute(planstartprimaryR);
		fftw_execute(planstartsecondaryL);
		fftw_execute(planstartsecondaryR);
	}
	else
	{
		fftw_execute(planprimarycurrL);
		fftw_execute(planprimarycurrR);
		fftw_execute(plansecondarycurrL);
		fftw_execute(plansecondarycurrR);
	}

	//destroy plans
	if (start == TRUE)
	{
		fftw_destroy_plan(planprimarycurrL);
		fftw_destroy_plan(planprimarycurrR);
		fftw_destroy_plan(plansecondarycurrL);
		fftw_destroy_plan(plansecondarycurrR);
		fftw_destroy_plan(planstartprimaryL);
		fftw_destroy_plan(planstartprimaryR);
		fftw_destroy_plan(planstartsecondaryL);
		fftw_destroy_plan(planstartsecondaryR);
	}
	else
	{
		fftw_destroy_plan(planprimarycurrL);
		fftw_destroy_plan(planprimarycurrR);
		fftw_destroy_plan(plansecondarycurrL);
		fftw_destroy_plan(plansecondarycurrR);
	}

	//AUTOEQ ALGORITHM

	//create inverse plans
	if (start == TRUE)
	{
		outplanprimarycurrL = fftw_plan_dft_1d(FRAMES_PER_BUFFER, outcbuffprimaryL, savecbuffprimaryL, FFTW_BACKWARD, flags);
		outplanprimarycurrR = fftw_plan_dft_1d(FRAMES_PER_BUFFER, outcbuffprimaryR, savecbuffprimaryR, FFTW_BACKWARD, flags);
		outplansecondarycurrL = fftw_plan_dft_1d(FRAMES_PER_BUFFER, outcbuffsecondaryL, savecbuffsecondaryL, FFTW_BACKWARD, flags);
		outplansecondarycurrR = fftw_plan_dft_1d(FRAMES_PER_BUFFER, outcbuffsecondaryR, savecbuffsecondaryR,  FFTW_BACKWARD, flags);
		outplanstartprimaryL = fftw_plan_dft_1d(FRAMES_PER_BUFFER/2, outcstartprimaryL, cstartprimaryL, FFTW_BACKWARD, flags);
		outplanstartprimaryR = fftw_plan_dft_1d(FRAMES_PER_BUFFER/2, outcstartprimaryR, cstartprimaryR, FFTW_BACKWARD, flags);
		outplanstartsecondaryL = fftw_plan_dft_1d(FRAMES_PER_BUFFER/2, outcstartsecondaryL, cstartsecondaryL, FFTW_BACKWARD, flags);
		outplanstartsecondaryR = fftw_plan_dft_1d(FRAMES_PER_BUFFER/2, outcstartsecondaryR, cstartsecondaryR, FFTW_BACKWARD, flags);
	}
	else
	{
		outplanprimarycurrL = fftw_plan_dft_1d(FRAMES_PER_BUFFER, outcbuffprimaryL, savecbuffprimaryL, FFTW_BACKWARD, flags);
		outplanprimarycurrR = fftw_plan_dft_1d(FRAMES_PER_BUFFER, outcbuffprimaryR, savecbuffprimaryR, FFTW_BACKWARD, flags);
		outplansecondarycurrL = fftw_plan_dft_1d(FRAMES_PER_BUFFER, outcbuffsecondaryL, savecbuffsecondaryL, FFTW_BACKWARD, flags);
		outplansecondarycurrR = fftw_plan_dft_1d(FRAMES_PER_BUFFER, outcbuffsecondaryR, savecbuffsecondaryR,  FFTW_BACKWARD, flags);
	}

	//execute inverse plans
	if (start == TRUE)
	{
		fftw_execute(outplanprimarycurrL);
		fftw_execute(outplanprimarycurrR);
		fftw_execute(outplansecondarycurrL);
		fftw_execute(outplansecondarycurrR);
		fftw_execute(outplanstartprimaryL);
		fftw_execute(outplanstartprimaryR);
		fftw_execute(outplanstartsecondaryL);
		fftw_execute(outplanstartsecondaryR);
	}
	else
	{		
		fftw_execute(outplanprimarycurrL);
		fftw_execute(outplanprimarycurrR);
		fftw_execute(outplansecondarycurrL);
		fftw_execute(outplansecondarycurrR);
	}

	//destroy plans
	if (start == TRUE)
	{
		fftw_destroy_plan(outplanprimarycurrL);
		fftw_destroy_plan(outplanprimarycurrR);
		fftw_destroy_plan(outplansecondarycurrL);
		fftw_destroy_plan(outplansecondarycurrR);
		fftw_destroy_plan(outplanstartprimaryL);
		fftw_destroy_plan(outplanstartprimaryR);
		fftw_destroy_plan(outplanstartsecondaryL);
		fftw_destroy_plan(outplanstartsecondaryR);
	}
	else
	{
		fftw_destroy_plan(outplanprimarycurrL);
		fftw_destroy_plan(outplanprimarycurrR);
		fftw_destroy_plan(outplansecondarycurrL);
		fftw_destroy_plan(outplansecondarycurrR);
	}

	//window and scale and output real output
	if (start == TRUE)
	{

		for (i = 0; i < FRAMES_PER_BUFFER/2; i++)
		{
			out[i*2] = (coswin[i] * cstartprimaryL[i][0] / (FRAMES_PER_BUFFER/2.0f)) + (coswin[i] * cstartsecondaryL[i][0] / (FRAMES_PER_BUFFER/2.0f)) + (win[i] * savecbuffprimaryL[i][0] / (float)FRAMES_PER_BUFFER) + (win[i] * savecbuffsecondaryL[i][0] / (float)FRAMES_PER_BUFFER);
			out[i*2+1] = (coswin[i] * cstartprimaryR[i][0] / (FRAMES_PER_BUFFER/2.0f)) + (coswin[i] * cstartsecondaryR[i][0] / (FRAMES_PER_BUFFER/2.0f)) + (win[i] * savecbuffprimaryR[i][0] / (float)FRAMES_PER_BUFFER) + (win[i] * savecbuffsecondaryR[i][0] / (float)FRAMES_PER_BUFFER);
		}
	}
	else
	{
		for (i = 0; i < FRAMES_PER_BUFFER/2; i++)
		{
			out[i*2] = (coswin[i] * primaryprevL[i+(FRAMES_PER_BUFFER/2)][0] / (float)FRAMES_PER_BUFFER) + (coswin[i] * secondaryprevL[i+(FRAMES_PER_BUFFER/2)][0] / (float)FRAMES_PER_BUFFER) + (win[i] * savecbuffprimaryL[i][0] / (float)FRAMES_PER_BUFFER) + (win[i] * savecbuffsecondaryL[i][0] / (float)FRAMES_PER_BUFFER);
			out[i*2+1] = (coswin[i] * primaryprevR[i+(FRAMES_PER_BUFFER/2)][0] / (float)FRAMES_PER_BUFFER) + (coswin[i] * secondaryprevR[i+(FRAMES_PER_BUFFER/2)][0] / (float)FRAMES_PER_BUFFER) + (win[i] * savecbuffprimaryR[i][0] / (float)FRAMES_PER_BUFFER) + (win[i] * savecbuffsecondaryR[i][0] / (float)FRAMES_PER_BUFFER);
		}
	}

	//free buffers
	if (start == TRUE)
	{
		fftw_free(cbuffprimaryL);
		fftw_free(cbuffprimaryR);
		fftw_free(outcbuffprimaryL);
		fftw_free(outcbuffprimaryR);
		fftw_free(cbuffsecondaryL);
		fftw_free(cbuffsecondaryR);
		fftw_free(outcbuffsecondaryL);
		fftw_free(outcbuffsecondaryR);
		fftw_free(realprimarycurr);
		fftw_free(realprimarycurrL);
		fftw_free(realprimarycurrR);
		fftw_free(realsecondarycurr);
		fftw_free(realsecondarycurrL);
		fftw_free(realsecondarycurrR);

		fftw_free(instartprimary);
		fftw_free(instartprimaryL);
		fftw_free(instartprimaryR);
		fftw_free(instartsecondary);
		fftw_free(instartsecondaryL);
		fftw_free(instartsecondaryR);
		fftw_free(cstartprimaryL);
		fftw_free(cstartprimaryR);
		fftw_free(outcstartprimaryL);
		fftw_free(outcstartprimaryR);
		fftw_free(cstartsecondaryL);
		fftw_free(cstartsecondaryR);
		fftw_free(outcstartsecondaryL);
		fftw_free(outcstartsecondaryR);
		start = FALSE;
	}
	else
	{
		fftw_free(cbuffprimaryL);
		fftw_free(cbuffprimaryR);
		fftw_free(outcbuffprimaryL);
		fftw_free(outcbuffprimaryR);
		fftw_free(cbuffsecondaryL);
		fftw_free(cbuffsecondaryR);
		fftw_free(outcbuffsecondaryL);
		fftw_free(outcbuffsecondaryR);
		fftw_free(realprimarycurr);
		fftw_free(realprimarycurrL);
		fftw_free(realprimarycurrR);
		fftw_free(realsecondarycurr);
		fftw_free(realsecondarycurrL);
		fftw_free(realsecondarycurrR);
	}

	//free buffers
	if (primaryprevL != NULL)
	{
		fftw_free(primaryprevL);
	}
	if (primaryprevR != NULL)
	{
		fftw_free(primaryprevR);
	}
	if (secondaryprevL != NULL)
	{
		fftw_free(secondaryprevL);
	}
	if (secondaryprevR != NULL)
	{
		fftw_free(secondaryprevR);
	}

	//point to saved buffers
	primaryprevL = savecbuffprimaryL;
	primaryprevR = savecbuffprimaryR;
	secondaryprevL = savecbuffsecondaryL;
	secondaryprevR = savecbuffsecondaryR;

	}
	return 0;
}
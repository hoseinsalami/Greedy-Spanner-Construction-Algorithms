#include <QtOpenGL>
#include <GL/glu.h>

#include "text.h"

struct vertex
{
	float x, y;      //Vertex
	float s, t;      //Texcoord
	vertex() {}
	vertex(float x, float y, float s, float t) : x(x), y(y), s(s), t(t) {}
};

static const int characters_length=512;
static vertex characters[characters_length];

static void logerr(const char * file, int line, const char * code){
	int err=glGetError();
	if (err!=GL_NO_ERROR)
		qDebug("GL Error at %s:%d '%s': %s", file, line, code, gluErrorString(err));
}
#define logerr(s) s;logerr(__FILE__,__LINE__,#s);
void a(){}

void write_text(float x, float y, QString line) {
	int sx=x;
	int len=line.length();
	int p=0;
	
	glPushMatrix();
	glLoadIdentity();

	//logerr(a());
	//logerr(glBindBuffer(GL_ARRAY_BUFFER, buffer));
	glEnableClientState(GL_VERTEX_ARRAY);
	glEnableClientState(GL_TEXTURE_COORD_ARRAY);
	logerr(glVertexPointer(2, GL_FLOAT, sizeof(vertex), &characters[0].x));
	logerr(glTexCoordPointer(2, GL_FLOAT, sizeof(vertex), &characters[0].s));

	for (int i=0; i<len; i++) {
		char c = line[i].toAscii();
		switch(c) {
			case ' ':  // Handle space separately
				x+=6;
				break;
			case '\n':  // Support newline character
				y+=16;
				x=sx;
				break;
			default:
				float cx = (c % 16) / 16.0f;
				float cy = (c / 16 - 2) / 8.0f;
		
				characters[p++]=vertex(x   , y+16, cx        , 1-cy-1/8.0f);
				characters[p++]=vertex(x+16, y+16, cx+1/16.0f, 1-cy-1/8.0f);
				characters[p++]=vertex(x+16, y   , cx+1/16.0f, 1-cy       );
				characters[p++]=vertex(x   , y   , cx        , 1-cy       );

				x+=8;
		}
		if (p>=characters_length) {
			//logerr(glBufferData(GL_ARRAY_BUFFER, sizeof(characters), characters, GL_STREAM_DRAW));
			logerr(glDrawArrays(GL_QUADS, 0, characters_length));
			p=0;
		}
	}
	if (p>0) {
		//logerr(glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(vertex)*p, characters));
		logerr(glDrawArrays(GL_QUADS, 0, p));
	}
	glDisableClientState(GL_TEXTURE_COORD_ARRAY);
	glDisableClientState(GL_VERTEX_ARRAY);
	
	glPopMatrix();
}


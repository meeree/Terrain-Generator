#include "/home/jhazelden/Cpp/OpenGL/loadinshader.cpp"
#define GLM_FORCE_RADIANS
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/ext.hpp>

#include <map>

#define SZ(n) ((n-1)*(n-1)*2*(3+2))

typedef std::vector<std::vector<GLfloat>> gridType;

struct Vertex 
{
    glm::vec3 pos; 
    GLfloat normal [3];
};

static struct 
{
	GLuint vao;
    GLint width, height;
    GLuint shaderProgram;
    GLuint pMatLoc, vMatLoc, mMatLoc;
    GLFWwindow* window;
    glm::vec3 camPos = glm::vec3(0,7,23);
    glm::vec3 camDir = glm::vec3(0,0,0);

    gridType grid;

    GLfloat corners [4] {0.0, 0.0, 0.0, 0.0};
    std::map<int, unsigned> keyMap {{GLFW_KEY_V,0}, {GLFW_KEY_B,1}, {GLFW_KEY_N,2}, {GLFW_KEY_M,3}};
    bool redraw;
    double randMax = 0.0, roughness = 0.5;
    double dist = 0.025;
    unsigned n = 1025;
    double lowest = 0.0, highest = 0.0;
    double weight = 0.0;

    std::vector<Vertex> verts;
    unsigned coloring = 0;
    unsigned colorCnt = 5;
} attrs;

double modfunc (double const& p, double const& q)
{
    int inc = std::signbit(p) == true ? -1 : 1;

    int i = 0;
    while (i*q <= fabs(p)) {++i;}
    return inc*i;
}

void updateVertex (Vertex& v, GLfloat const& p1, GLfloat const& p2, GLfloat const& p3) 
{
    v.pos[0] = attrs.dist*p1; 
    v.pos[1] = p2; 
//    v.pos[1] = 0.1*modfunc(p2, 0.1);
    v.pos[2] = attrs.dist*p3;
}

void constructNormals (Vertex& v1, Vertex& v2, Vertex& v3)
{
    glm::vec3 p1 {v1.pos[0], v1.pos[1], v1.pos[2]};
    glm::vec3 p2 {v2.pos[0], v2.pos[1], v2.pos[2]};
    glm::vec3 p3 {v3.pos[0], v3.pos[1], v3.pos[2]};
    glm::vec3 nm {glm::normalize(glm::cross(p2-p1, p3-p1))};

    v1.normal[0] = nm[0]; v1.normal[1] = nm[1]; v1.normal[2] = nm[2];
    v2.normal[0] = nm[0]; v2.normal[1] = nm[1]; v2.normal[2] = nm[2];
    v3.normal[0] = nm[0]; v3.normal[1] = nm[1]; v3.normal[2] = nm[2];
}

void setShaders (char const* vertLoc, char const* fragLoc)
{
	attrs.shaderProgram = glCreateProgram();

    auto vertShader = loadInShader(vertLoc, GL_VERTEX_SHADER);
    auto fragShader = loadInShader(fragLoc, GL_FRAGMENT_SHADER);

    glAttachShader(attrs.shaderProgram, vertShader);
    glAttachShader(attrs.shaderProgram, fragShader);

    glDeleteShader(vertShader);
    glDeleteShader(fragShader);

    glLinkProgram(attrs.shaderProgram);
}

void startup (GLfloat const& width, GLfloat const& height, char const* vertLoc, char const* fragLoc, const char* title="Untitled Window")
{
	attrs.width = width; 
    attrs.height = height;

    if(!glfwInit()) {
        std::cerr<<"failed to initialize glfw"<<std::endl;
        exit(EXIT_SUCCESS);
    }
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 5);
    attrs.window = glfwCreateWindow(attrs.width, attrs.height, title, NULL, NULL);

    if(!attrs.window) {
        std::cerr<<"failed to initialize window"<<std::endl;
        exit(EXIT_SUCCESS);
    }
    glfwMakeContextCurrent(attrs.window);

    glewExperimental = GL_TRUE;
    if(glewInit() != 0) {
        std::cerr<<"failed to initialize glew"<<std::endl;
        exit(EXIT_SUCCESS);
    }

    setShaders(vertLoc, fragLoc);
    glUseProgram(attrs.shaderProgram);
}

void shutdown () 
{
	glDeleteVertexArrays(1, &attrs.vao);
    glDeleteProgram(attrs.shaderProgram);
}

void setTitle () 
{
    glfwSetWindowTitle(attrs.window, 
            ("rand max: "+std::to_string(attrs.randMax)+", roughness: "+std::to_string(attrs.roughness)+", weight: "+std::to_string(attrs.weight)
             +", crnr. 0:  "+std::to_string(attrs.corners[0])+", crnr. 1: "+std::to_string(attrs.corners[1])+", crnr. 2: "+std::to_string(attrs.corners[2])+", crnr. 3: "+std::to_string(attrs.corners[3])).c_str());
}

void resize ()
{
    attrs.verts.resize(0);
    attrs.grid.resize(0);
    for (auto& v: attrs.grid) 
    {
	    v.resize(0);
    }

    attrs.verts.resize(SZ(attrs.n));
    attrs.grid.resize(attrs.n);
    for (auto& v: attrs.grid) 
    {
	    v.resize(attrs.n);
    }
}

bool isEmpty (std::fstream& fl)
{
    return fl.peek() == std::ifstream::traits_type::eof();
}

void saveState (std::fstream& fl)
{
    fl<<std::to_string(attrs.randMax)<<";\n"
    <<std::to_string(attrs.roughness)<<";\n"
    <<std::to_string(attrs.weight)<<";\n"
    <<std::to_string(attrs.corners[0])<<","
    <<std::to_string(attrs.corners[1])<<","
    <<std::to_string(attrs.corners[2])<<","
    <<std::to_string(attrs.corners[3])<<";\n";
}

void getState (std::fstream& fl)
{
    if (isEmpty(fl)) 
    {
        saveState(fl);
        return;
    }
    fl>>attrs.randMax; fl.ignore();
    fl>>attrs.roughness; fl.ignore();
    fl>>attrs.weight; fl.ignore();
    fl>>attrs.corners[0]; fl.ignore();
    fl>>attrs.corners[1]; fl.ignore();
    fl>>attrs.corners[2]; fl.ignore();
    fl>>attrs.corners[3];
}

void keyCallback(GLFWwindow*, int key, int, int action, int)
{   
    if (key == GLFW_KEY_S && (action == GLFW_REPEAT || action == GLFW_PRESS))
        attrs.camPos[0] -= 0.1f;
    if (key == GLFW_KEY_W && (action == GLFW_REPEAT || action == GLFW_PRESS))
        attrs.camPos[0] += 0.1f;
    if (key == GLFW_KEY_Q && (action == GLFW_REPEAT || action == GLFW_PRESS))
        attrs.camPos[1] -= 0.1f;
    if (key == GLFW_KEY_E && (action == GLFW_REPEAT || action == GLFW_PRESS))
        attrs.camPos[1] += 0.1f;
    if (key == GLFW_KEY_A && (action == GLFW_REPEAT || action == GLFW_PRESS))
        attrs.camPos[2] -= 0.1f;
    if (key == GLFW_KEY_D && (action == GLFW_REPEAT || action == GLFW_PRESS))
        attrs.camPos[2] += 0.1f;
    if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS)
    {
        std::fstream fl;
        fl.open("/home/jhazelden/Cpp/OpenGL/diamondSquare/fl", std::ios::out | std::ios::trunc);
        saveState(fl);
        fl.close();

	    glfwTerminate(); 
        exit(EXIT_SUCCESS);
    }
    if (key == GLFW_KEY_R && action == GLFW_PRESS)
    {
        attrs.redraw = true;
    }
    if (key == GLFW_KEY_K && action == GLFW_PRESS)
    {
        attrs.coloring = (attrs.coloring+1)%attrs.colorCnt;
    }
    if (glfwGetKey(attrs.window, GLFW_KEY_J) == GLFW_PRESS)
        attrs.camDir += glm::vec3(0.1f, 0.0f, 0.0f);
    if (glfwGetKey(attrs.window, GLFW_KEY_L) == GLFW_PRESS)
        attrs.camDir += glm::vec3(-0.1f, 0.0f, 0.0f);
    if (glfwGetKey(attrs.window, GLFW_KEY_UP) == GLFW_PRESS) 
    {
        for (auto const& keyPair: attrs.keyMap)
        {
            if (glfwGetKey(attrs.window, keyPair.first) == GLFW_PRESS)
            {
                attrs.corners[keyPair.second] += 10*attrs.dist;
            }
        }
        if (glfwGetKey(attrs.window, GLFW_KEY_Z) == GLFW_PRESS)
        {
            attrs.randMax += attrs.dist;
        }
        if (glfwGetKey(attrs.window, GLFW_KEY_X) == GLFW_PRESS)
        {
            attrs.roughness = fmod(attrs.roughness+0.05, 1.0);
        }
        if (glfwGetKey(attrs.window, GLFW_KEY_C) == GLFW_PRESS)
        {
            for (auto& corner: attrs.corners)
            {
                corner /= attrs.dist;
            }
            attrs.randMax /= attrs.dist;

            attrs.dist *= 10./9;

            for (auto& corner: attrs.corners)
            {
                corner *= attrs.dist;
            }
            attrs.randMax *= attrs.dist;
        }
        if (glfwGetKey(attrs.window, GLFW_KEY_L) == GLFW_PRESS)
        {
            attrs.n = (attrs.n-1)*2+1;
        }
        if (glfwGetKey(attrs.window, GLFW_KEY_H) == GLFW_PRESS)
        {
            attrs.weight += 0.1;
        }
        setTitle();
    }
    else if (glfwGetKey(attrs.window, GLFW_KEY_DOWN) == GLFW_PRESS) 
    {
        for (auto const& keyPair: attrs.keyMap)
        {
            if (glfwGetKey(attrs.window, keyPair.first) == GLFW_PRESS)
            {
                attrs.corners[keyPair.second] -= 10*attrs.dist;
            }
        }
        if (glfwGetKey(attrs.window, GLFW_KEY_Z) == GLFW_PRESS)
        {
            attrs.randMax -= attrs.dist;
        }
        if (glfwGetKey(attrs.window, GLFW_KEY_X) == GLFW_PRESS)
        {
            attrs.roughness = fmod(attrs.roughness-0.05, 1.0);
        }
        if (glfwGetKey(attrs.window, GLFW_KEY_C) == GLFW_PRESS)
        {
            for (auto& corner: attrs.corners)
            {
                corner /= attrs.dist;
            }
            attrs.randMax /= attrs.dist;

            attrs.dist *= 9./10;

            for (auto& corner: attrs.corners)
            {
                corner *= attrs.dist;
            }
            attrs.randMax *= attrs.dist;
        }
        if (glfwGetKey(attrs.window, GLFW_KEY_L) == GLFW_PRESS)
        {
            attrs.n = (attrs.n-1)/2+1;
        }
        if (glfwGetKey(attrs.window, GLFW_KEY_H) == GLFW_PRESS)
        {
            attrs.weight -= 0.1;
        }
        setTitle();
    }
}

void redraw ();

void render (double const& t, double const&) 
{   
    if (attrs.redraw)
    {
        resize();
        redraw();
    }
    
    GLfloat const color [4] {0.0f, 0.0f, 0.0f, 1.0f};
    glClearBufferfv(GL_COLOR, 0.0f, color);
    glClear(GL_DEPTH_BUFFER_BIT);
    glUniform1f(2, t);
    glUniform1ui(3, attrs.coloring);
    glUniform1f(4, attrs.lowest);
    glUniform1f(5, attrs.highest);
    glUniform1f(6, attrs.dist);

//    auto vMat = glm::lookAt(attrs.camPos, attrs.camDir, glm::normalize(glm::cross(attrs.camPos, attrs.camDir))); 
    auto vMat = glm::lookAt(attrs.camPos, attrs.camDir, glm::vec3(0.0, 1.0, 0.0)); 
    glUniformMatrix4fv(attrs.vMatLoc, 1, GL_FALSE, glm::value_ptr(vMat));
//    auto angle = glm::radians((GLfloat)t*3.0f);
//    auto mMat = glm::rotate(angle, glm::vec3(0.0f, 1.0f, 0.0f));
    glDrawArrays(GL_TRIANGLES, 0, 3*SZ(attrs.n)/5);
    glDrawArrays(GL_LINES, 3*SZ(attrs.n)/5, 2*SZ(attrs.n)/5);

    glfwPollEvents();
    glfwSwapBuffers(attrs.window);
}

void loop (double t=0.0, double const& dt=0.1) 
{
	while (true)
    {
	    render(t, dt);
        t += dt;
    }
}

bool gridCheck (unsigned const& n, gridType const& grid) 
{
	return grid.size() == n && std::all_of(grid.begin(), grid.end(), [&](std::vector<GLfloat> const& v){return v.size() == n;});
}

void tesselate (unsigned const& n, gridType const& grid, std::vector<Vertex>& verts)
{
	if (!gridCheck(n, grid)) 
    {
	    std::cerr<<"grid not"<<n<<"X"<<n<<"."<<std::endl;
        exit(EXIT_SUCCESS);
    }    
//    if (verts.size() != (n-1)*(n-1)*2*3)
//    {
//	    std::cerr<<"verts size incorrect"<<std::endl;
//        exit(EXIT_SUCCESS);
//    }

    int k = 0;
    for (int i = 0; i < (int)n-1; ++i)
    {
        for (int j = 0; j < (int)n-1; ++j, k += 6) 
        {
	        int p {((int)n-1)/2};

            updateVertex(verts[k], j-p, grid[i+1][j], i+1-p);
            updateVertex(verts[k+1], j+1-p, grid[i][j+1], i-p);
            updateVertex(verts[k+2], j-p, grid[i][j], i-p);
            constructNormals(verts[k], verts[k+1], verts[k+2]);

            updateVertex(verts[k+3], j+1-p, grid[i][j+1], i-p);
            updateVertex(verts[k+4], j-p, grid[i+1][j], i+1-p);
            updateVertex(verts[k+5], j+1-p, grid[i+1][j+1], i+1-p);
            constructNormals(verts[k+3], verts[k+4], verts[k+5]);
        }
    }
}

void maxCheck (GLfloat& v)
{   
    if (v >= attrs.highest*0.7)
    {
        v = attrs.highest*0.7+attrs.dist*5*(v/attrs.highest);
    }
}

void minCheck (GLfloat& v)
{   
    float k = (v-attrs.lowest)/(attrs.highest-attrs.lowest);
    if (k < 0.3)
    {
        v = 0.3*(attrs.highest-attrs.lowest)+attrs.lowest;
    }
}

void diamondStep (unsigned const& n, gridType& grid, std::pair<unsigned, unsigned>* iter, unsigned const& k, std::mt19937& gen, std::uniform_real_distribution<>& dis, bool const& grow)
{
	for (unsigned i = 0; i < ((n-1)/k)*((n-1)/k); ++i)
    {
	    auto p = *iter++;

        if (grow)
        {
            grid[p.first][p.second] +=
            (grid[p.first-k/2][p.second-k/2]
            +grid[p.first-k/2][p.second+k/2]
            +grid[p.first+k/2][p.second-k/2]
            +grid[p.first+k/2][p.second+k/2])/4
            +dis(gen);
        }
        else
        {
            grid[p.first][p.second] =
            (grid[p.first-k/2][p.second-k/2]
            +grid[p.first-k/2][p.second+k/2]
            +grid[p.first+k/2][p.second-k/2]
            +grid[p.first+k/2][p.second+k/2])/4
            +dis(gen);
        }
    }
}

void squareStep (unsigned const& n, gridType& grid, std::pair<unsigned, unsigned>* iter, unsigned const& k, std::mt19937& gen, std::uniform_real_distribution<>& dis, bool const& grow)
{
	for (unsigned i = 0; i < ((n-1)/k)*((n-1)/k); ++i)
    {
	    auto p = *iter++;
        if (p.first == n-k/2-1 || p.second == n-k/2-1)
        {
            continue;
        }
        auto overlapx = (p.first == n-k/2-1 ? grid[k/2][p.second] : grid[p.first+k][p.second]);
        auto overlapz = (p.second == n-k/2-1 ? grid[p.first][k/2] : grid[p.first][p.second+k]);

        if (grow)
        {
            grid[p.first+k/2][p.second] +=
            (grid[p.first][p.second]
            +grid[p.first+k/2][p.second-k/2]
            +grid[p.first-k/2][p.second-k/2]
            +overlapz)/4
            +dis(gen);

            grid[p.first][p.second+k/2] +=
            (grid[p.first][p.second]
            +grid[p.first+k/2][p.second+k/2]
            +grid[p.first-k/2][p.second-k/2]
            +overlapx)/4
            +dis(gen);
        }
        else
        {
            grid[p.first+k/2][p.second] =
            (grid[p.first][p.second]
            +grid[p.first+k/2][p.second-k/2]
            +grid[p.first-k/2][p.second-k/2]
            +overlapz)/4
            +dis(gen);

            grid[p.first][p.second+k/2] =
            (grid[p.first][p.second]
            +grid[p.first+k/2][p.second+k/2]
            +grid[p.first-k/2][p.second-k/2]
            +overlapx)/4
            +dis(gen);
        }

    }
}

void diamondSquare (unsigned const& n, gridType& grid, double min, double max, double decreaseFactor, bool const& grow)
{
	double realPart;
    if (modf(log2(n-1), &realPart) > 0.01) 
    {
	    std::cerr<<"n for diamond square not in format 2^k+1"<<std::endl;
        exit(EXIT_SUCCESS);
    } 
    
    std::random_device rd;
    std::mt19937 gen {rd()};

    for (unsigned p=log2(n-1); p > 0; --p)
    {
	    unsigned k = pow(2,p);
        std::pair<unsigned, unsigned> points [((n-1)/k)*((n-1)/k)];
        for (int i = 0; i < ((int)n-1)/(int)k; ++i)
            for (int j = 0; j < ((int)n-1)/(int)k; ++j)
            {
	            unsigned z = i*k+k/2;
                unsigned x = j*k+k/2;
                points[i*((n-1)/k)+j] = {z,x};
            }
        std::uniform_real_distribution<> dis (min, max);
        diamondStep(n, grid, points, k, gen, dis, grow);
        squareStep(n, grid, points, k, gen, dis, grow);
        min *= decreaseFactor;
        max *= decreaseFactor;
    }
}

void smooth (unsigned const& n, gridType& grid, double const& weight)
{   
    if (!gridCheck(n, grid)) 
    {
	    std::cerr<<"grid not"<<n<<"X"<<n<<"."<<std::endl;
        exit(EXIT_SUCCESS);
    }    

    gridType newGrid (n);
    for (auto& v: newGrid)
    {
        v.resize(n);
        std::fill(v.begin(), v.end(), 0.0);
    }

    for (int i = 1; i < (int)n-1; ++i)
    {
	    for (int j = 1; j < (int)n-1; ++j)
        {
#define SMOOTH9

#ifndef SMOOTH4 
            newGrid[i][j] = 
                (weight*grid[i][j]
                +grid[i-1][j]+grid[i+1][j]+grid[i][j-1]+grid[i][j+1])/(4*weight);
#endif //SMOOTH4

#ifndef SMOOTH9
            weightp = (attrs.highest-grid[i][j])/attrs.highest;
            newGrid[i][j] = 
                (weightp*grid[i][j]
                +grid[i-1][j-1]+grid[i-1][j]+grid[i-1][j+1]
                +grid[i][j-1]+grid[i][j]+grid[i][j+1]
                +grid[i+1][j-1]+grid[i+1][j]+grid[i+1][j+1])/(9+weightp);
#endif //SMOOTH9
        }
    }

    std::copy(newGrid.begin(), newGrid.end(), grid.begin());
}

void setMinMax (gridType const& grid)
{
    attrs.highest = 0.0;
    attrs.lowest = 0.0;
    for (auto const& vec: grid)
    {
        for (auto const& v: vec )
        {
            if (v > attrs.highest)
            {
                attrs.highest = v;
            }
            else if (v < attrs.lowest)
            {
                attrs.lowest = v; 
            }
        }
    }
}

void addGrass ()
{
    unsigned j = 3*SZ(attrs.n)/5;
	for (unsigned i = 0; i < 3*SZ(attrs.n)/5; i += 3, j += 2)
	{
		Vertex const& v1 {attrs.verts[i]};
		Vertex const& v2 {attrs.verts[i+1]};
		Vertex const& v3 {attrs.verts[i+2]};
        auto avgHeight = (v1.pos[1]+v2.pos[1]+v3.pos[1])/3;
        float k = (avgHeight-attrs.lowest)/(attrs.highest-attrs.lowest);
#define THRESHOLD 0.001
        if (k <= 0.3+THRESHOLD)
        {
            continue;
        }
		
		glm::vec3 nm {glm::vec3(v1.normal[0], v1.normal[1], v1.normal[2])};
		double theta {glm::acos(glm::dot(nm, glm::vec3(0.0, 1.0, 0.0)))};

		if (theta < 3.1415/7)
		{
            Vertex l1 {{(v1.pos[0]+v2.pos[0]+v3.pos[0])/3, (v1.pos[1]+v2.pos[1]+v3.pos[1])/3, (v1.pos[2]+v2.pos[2]+v3.pos[2])/3}, {2,0,0}};
            GLfloat scl = (3.14/2-theta)*attrs.dist;
            Vertex l2 {{l1.pos[0]+scl*nm[0], l1.pos[1]+scl*nm[1], l1.pos[2]+scl*nm[2]}, {2,0,0}};
			attrs.verts[j] = l1;
			attrs.verts[j+1] = l2;
		}
	}
}

void redraw ()
{
    setTitle();
    attrs.grid[attrs.n-1][attrs.n-1] = attrs.corners[0];
    attrs.grid[0][attrs.n-1] = attrs.corners[1];
    attrs.grid[attrs.n-1][0] = attrs.corners[2];
    attrs.grid[0][0] = attrs.corners[3];

    diamondSquare(attrs.n, attrs.grid, -attrs.randMax, attrs.randMax, attrs.roughness, false);
    smooth(attrs.n, attrs.grid, attrs.weight);
    setMinMax(attrs.grid);
    for (auto& vec: attrs.grid)
    {
	    for (auto& v: vec)
        {
	        maxCheck(v);
            minCheck(v);
        }
    }

    tesselate(attrs.n, attrs.grid, attrs.verts);     
    addGrass();
    glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(Vertex)*attrs.verts.size(), attrs.verts.data());

    attrs.redraw = false;
}

int main ()
{
    {
        std::fstream fl;
        fl.open("/home/jhazelden/Cpp/OpenGL/diamondSquare/fl", std::ios::in);
        getState(fl);
        fl.close();

        resize();

        startup(1920, 1080, "/home/jhazelden/Cpp/OpenGL/diamondSquare/Shaders/vert.glsl", "/home/jhazelden/Cpp/OpenGL/diamondSquare/Shaders/frag.glsl", "Terrain Generation Practice"); 

        attrs.redraw = true;

        GLuint vao;
        glGenVertexArrays(1, &vao);
        glBindVertexArray(vao);

        GLuint buffer;
        glGenBuffers(1, &buffer);
        glBindBuffer(GL_ARRAY_BUFFER, buffer);
        glBufferData(GL_ARRAY_BUFFER, sizeof(Vertex)*attrs.verts.size(), attrs.verts.data(), GL_DYNAMIC_DRAW);
        
        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex), reinterpret_cast<void*>(offsetof(Vertex, pos)));
        glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex), reinterpret_cast<void*>(offsetof(Vertex, normal)));
        glEnableVertexAttribArray(0);
        glEnableVertexAttribArray(1);

        attrs.pMatLoc = glGetUniformLocation(attrs.shaderProgram, "pMat");
        attrs.vMatLoc = glGetUniformLocation(attrs.shaderProgram, "vMat");
        attrs.mMatLoc = glGetUniformLocation(attrs.shaderProgram, "mMat");

        // Init mat4s
        auto pMat = glm::perspective(glm::radians(45.0f),(GLfloat)attrs.width/attrs.height, 0.1f, 100.0f); 
        auto mMat = glm::mat4(1.0f);
        glUniformMatrix4fv(attrs.pMatLoc, 1, GL_FALSE, glm::value_ptr(pMat));
        glUniformMatrix4fv(attrs.mMatLoc, 1, GL_FALSE, glm::value_ptr(mMat));
        glUniformMatrix4fv(attrs.vMatLoc, 1, GL_FALSE, glm::value_ptr(glm::mat4(1.0f)));

        glfwSetKeyCallback(attrs.window, keyCallback);

        glLineWidth(10);
        glEnable (GL_BLEND); glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        glEnable(GL_DEPTH_TEST); glDepthFunc(GL_LESS);
        glEnable(GL_CULL_FACE); glCullFace(GL_BACK);
    }

    loop();
}

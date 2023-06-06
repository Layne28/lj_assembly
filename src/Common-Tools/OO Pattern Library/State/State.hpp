class state {
public:
    state();
    void bind();
    void unbind();
    double force(receptor*);


};

class bound : state {

};

class unbound : state {

};

class receptor {

};



class RunAndTumbleState {
    virtual void update(RunAndTumbleParticle* );
};

class Run: RunAndTumbleState {

};
class Tumble: RunAndTumbleState {

};
class Particle;
class RunAndTumbleParticle: Particle { 
private:
    RunAndTumbleState _state;

    



};
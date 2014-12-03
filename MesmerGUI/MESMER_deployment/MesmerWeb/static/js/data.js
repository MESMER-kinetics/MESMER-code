var app = app || {};

app.constant = {
    reactionType: {
        'one2one': 'A -> B',
        'one2two': 'A -> B + C',
        'two2one': 'A + B -> C',
        'two2two': 'A + B -> C + D'
    },
    reactantType: {
        model: 'modelled',
        deficient: 'deficientReactant',
        excess: 'excessReactant',
        transition: 'transitionState',
        sink: 'sink'
    },
    moleculeType: {
        ER: 'excessReactant',
        DR: 'deficientReactant',
        SM: 'sink',
        MM: 'modelled',
        TS: 'transitionState',
        BG: 'bathGas'
    },
    moleculeParameter: {
        'excessReactant': ['MW', 'sigma', 'epsilon', 'deltaEDown', 'deltaEDownTExponent', 'referenceTemperature'],
        'deficientReactant': ['MW', 'sigma', 'epsilon', 'deltaEDown', 'deltaEDownTExponent', 'referenceTemperature'],
        'sink': ['MW', 'sigma', 'epsilon', 'rotConsts', 'symmetryNumber', 'vibFreqs', 'frequenciesScaleFactor', 'spinMultiplicity', 'deltaEDown', 'deltaEDownTExponent', 'referenceTemperature'],
        'modelled': [],
        'transitionState': ['MW', 'sigma', 'epsilon', 'deltaEDown', 'deltaEDownTExponent', 'referenceTemperature'],
        'bathGas': ['rotConsts', 'symmetryNumber', 'ZPE', 'vibFreqs', 'frequenciesScaleFactor', 'spinMultiplicity', 'deltaEDown', 'deltaEDownTExponent', 'referenceTemperature']
    }
};

/*Molecule Model*/
app.Molecule = Backbone.Model.extend({
    defaults: {
        type: app.constant.moleculeType.MM,
        id: '',
        ZPE: '',
        ZPE_unit: 'cm-1',
        rotConsts: '',
        rotConsts_unit: 'cm-1',
        vibFreqs: '',
        frequenciesScaleFactor: '',
        symmetryNumber: '',
        MW: '',
        spinMultiplicity: '',
        epsilon: '',
        sigma: '',
        deltaEDown: '',
        DOSCMethod: 'ClassicalRotors',
        deltaEDownTExponent: '',
        deltaEDown: '',
        referenceTemperature: '298'
    },

    urlRoot: '/',

    sync: function(method, model) {

    },

    fromJSON: function(data) {
        if (!data)
            return;
        this.set(data);
    },

    convert_unit: function(data) {
        // convert ZPE_unit to default(1/cm)
        switch (data['ZPE_unit']) {
            case 'kCal/mol':
                data['ZPE'] = data['ZPE'] * 349.757;
                break;
            case 'kJ/mol':
                data['ZPE'] = data['ZPE'] * 83.593;
                break;
            default:
                break;
        }
        data['ZPE_unit'] = 'cm-1';

        // convert rotConsts_unit to default(1/cm)
        switch (data['rotConsts_unit']) {
            case 'GHz':
                value = data['rotConsts'].split(/\s+/);
                for (var i = 0; i < value.length; i++) {
                    value[i] = parseFloat(value[i]);
                    value[i] *= 0.0334;
                }
                data['rotConsts'] = value.join(' ');
                break;
            default:
                break;
        }
        data['rotConsts_unit'] = 'cm-1'
    }
});

app.MoleculeList = Backbone.Collection.extend({
    model: app.Molecule,

    url: '/',

    seq: 0,

    initialize: function() {
        var local = this;
        $.get('/library/', function(response, status, xhr) {
            local.library = JSON.parse(response);
            local.trigger('library');
        })
    },

    sync: function(method, model) {

    },

    nextOrder: function() {
        return this.length + 1;
    },

    select: function(el) {
        this.selected = el;
        this.trigger("select", el);
        el.trigger("selected");
    },

    fromJSON: function(list) {
        /*if (list){
			for (i=0; i<list.length; i++){
				molecule = new app.Molecule(list[i]);
				this.add(molecule);
			}
		}*/
        if (list) {
            var new_models = [];
            for (i = 0; i < list.length; i++) {
                molecule = new app.Molecule(list[i]);
                new_models.push(molecule);
            }
            this.reset(new_models);
            //this.selected = null;
            if (new_models.length > 0)
                this.select(new_models[0]);
            else
                this.selected = null;
        }
    }
});

app.Reaction = Backbone.Model.extend({
    defaults: {
        id: '',
        type: app.constant.reactionType.one2one,
        R1Type: app.constant.reactantType.model,
        R1Ref: '',

        R2Type: app.constant.reactantType.model,
        R2Ref: '',

        P1Type: app.constant.reactantType.model,
        P1Ref: '',

        P2Type: app.constant.reactantType.model,
        P2Ref: '',

        TType: app.constant.reactantType.transition,
        TRef: '',

        MCRCMethod: 'SimpleRRKM',
        preExponential: '',
        preExponential_unit: 'cm3molecule-1s-1',
        activationEnergy: '',
        activationEnergy_unit: 'kJ/mol',
        nInfinity: '',
        excessReactantConc: ''
    },
    sync: function(method, model) {

    },
    urlRoot: '/',

    fromJSON: function(data) {
        if (!data)
            return;
        this.set(data);
    }
});

app.ReactionList = Backbone.Collection.extend({
    model: app.Reaction,

    seq: 0,

    nextOrder: function() {
        return this.length + 1;
    },
    sync: function(method, model) {

    },

    url: '/',

    select: function(el) {
        this.selected = el;
        this.trigger("select", el);
        el.trigger("selected");
    },

    fromJSON: function(list) {
        if (list) {
            var new_reactions = [];
            for (var i = 0; i < list.length; i++) {
                var reaction = new app.Reaction(list[i]);
                new_reactions.push(reaction);
            }
            this.reset(new_reactions);
            if (new_reactions)
                this.select(new_reactions[0]);
            else
                this.selected = null;
        }
    }
});

app.Control = Backbone.Model.extend({
    defaults: {
        eigenvalues: '',
        calcMethod: 'simpleCalc',
        options: {
            printCellDOS: false,
            printCellTransitionStateFlux: false,
            printReactionOperatorColumnSums: false,
            printGrainBoltzmann: false,
            printGrainDOS: false,
            printGrainkbE: false,
            printGrainkfE: false,
            printTSsos: false,
            printGrainedSpeciesProfile: false,
            printGrainTransitionStateFlux: false,
            printReactionOperatorSizematrix: false,
            printSpeciesProfile: false,
            printTunnelingCoefficients: false,
            printTunnellingCoefficients: false,
            printCrossingCoefficients: false,
            MaximumEvolutionTime: false,
            testDOS: false,
            testMicroRates: false,
            testRateConstants: false,
            useTheSameCellNumberForAllConditions: false,
            hideInactive: false
        }
    },
    sync: function(method, model) {

    },
    urlRoot: '/',

    fromJSON: function(data) {
        this.set(data);
    }
});

app.ModelParameters = Backbone.Model.extend({
    defaults: {
        grainSize: '',
        energyAboveTheTopHill: ''
    },
    sync: function(method, model) {

    },
    urlRoot: '/',

    fromJSON: function(data) {
        if (!data)
            return;
        this.set(data);
    }
});

app.PT = Backbone.Model.extend({
    defaults: {
        P: '',
        T: '',
        precision: 'double'
    },
    sync: function(method, model) {

    },
    urlRoot: '/',

    fromJSON: function(data) {
        if (!data)
            return;
        this.set(data);
    }
});

app.PTS = Backbone.Collection.extend({
    model: app.PT,
    sync: function(method, model) {

    },

    url: '/',

    fromJSON: function(list) {
        //if (list){
        //	for (i=0; i<list.length; i++){
        //		pt = new app.PT(list[i]);
        //		this.add(pt);
        //	}
        //}
        if (list) {
            var new_models = [];
            for (i = 0; i < list.length; i++) {
                pt = new app.PT(list[i]);
                new_models.push(pt);
            }
            this.reset(new_models);
        }
    }
});

app.Conditions = Backbone.Model.extend({
    defaults: {
        bathGas: '',
        units: 'PPCC',
        initialPopulation: '',
        initialPopulationRef: '',
        PTs: new app.PTS()
    },
    sync: function(method, model) {

    },
    urlRoot: '/',

    fromJSON: function(data) {
        if (data) {
            this.set('bathGas', data['bathGas'] || '');
            this.set('units', data['units'] || '');
            this.set('initialPopulation', data['initialPopulation'] || '');
            this.set('initialPopulationRef', data['initialPopulationRef'] || '');
            this.get('PTs').fromJSON(data['PTs']);
        }
    }
});

app.MesmerData = Backbone.Model.extend({
    defaults: {
        title: 'MESMER Title',
        moleculeList: new app.MoleculeList(),
        reactionList: new app.ReactionList(),
        control: new app.Control(),
        modelParameters: new app.ModelParameters(),
        conditions: new app.Conditions()
    },

    sync: function(method, model) {

    },

    urlRoot: '/',

    toJSON: function() {
        var result = {};

        result.title = this.get('title');
        var moleculeList = this.get('moleculeList').models;
        result.moleculeList = [];
        for (i = 0; i < moleculeList.length; i++) {
            var molecule = moleculeList[i].toJSON();
            //moleculeList[i].convert_unit(molecule);
            result.moleculeList.push(molecule);
        }

        var reactionList = this.get('reactionList').models;
        result.reactionList = [];
        for (i = 0; i < reactionList.length; i++)
            result.reactionList.push(reactionList[i].toJSON());

        var control = this.get('control');
        result.control = control.toJSON();

        var modelParameters = this.get('modelParameters');
        result.modelParameters = modelParameters.toJSON();

        var conditions = this.get('conditions');
        result.conditions = conditions.toJSON();

        return result;
    },

    fromJSON: function(data) {
        this.set('title', data['title'])
        this.get('moleculeList').fromJSON(data['moleculeList']);
        this.get('reactionList').fromJSON(data['reactionList']);
        this.get('control').fromJSON(data['control']);
        this.get('modelParameters').fromJSON(data['modelParameters']);
        this.get('conditions').fromJSON(data['conditions']);
    }
});
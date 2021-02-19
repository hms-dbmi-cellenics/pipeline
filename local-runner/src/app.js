const AWS = require('aws-sdk');
const fs = require('fs');
const Docker = require('dockerode');
const DockerEvents = require('docker-events');
const { pt } = require('prepend-transform');
const chalk = require('chalk');

// Enable & connect to local Inframock.
AWS.config.update({
  endpoint: 'http://localhost:4566',
  sslEnabled: false,
  s3ForcePathStyle: true,
});

const initStack = async () => {
  console.log('Loading CloudFormation for local container launcher...');
  const data = fs.readFileSync('cf-local-container-launcher.yaml', 'utf-8');

  console.log('Creating mock Lambda function on InfraMock...');
  const cf = new AWS.CloudFormation({
    region: 'eu-west-1',
  });

  try {
    await cf.deleteStack({
      StackName: 'local-container-launcher',
    }).promise();
  } catch (e) {
    console.log('No previous stack found on InfraMock.');
  }

  const { StackId } = await cf.createStack({
    StackName: 'local-container-launcher',
    Capabilities: ['CAPABILITY_IAM', 'CAPABILITY_NAMED_IAM'],
    TemplateBody: data,
  }).promise();

  console.log('Stack with ARN', StackId, 'successfully created.');
};

const generateRandomColor = (nameColorMap, name) => {
  // Generate random color for output and easier understanding.
  const randomHexValue = () => Math.floor(Math.random() * 255);
  const coloredName = chalk.rgb(randomHexValue(), randomHexValue(), randomHexValue())(name);

  // eslint-disable-next-line no-param-reassign
  nameColorMap[name] = coloredName;

  return coloredName;
};

const attachToExistingContainers = (docker, nameColorMap) => {
  docker.listContainers((err, containers) => {
    containers.forEach((info) => {
      const { Names: names, Id: id } = info;
      const name = names.filter((n) => n.includes('pipeline'))[0];

      if (!name) {
        return;
      }

      const coloredName = generateRandomColor(nameColorMap, name.replace('/', ''));
      console.log('Container with name', coloredName, 'already running, attaching log...');

      docker.getContainer(id).attach(
        { stream: true, stdout: true, stderr: true },
        (err, stream) => {
          stream.pipe(pt(`${coloredName} | `)).pipe(process.stdout);
        },
      );
    });
  });
};

const attachToNewContainers = (docker, emitter, nameColorMap) => {
  emitter.on('start', (message) => {
    const { id, Actor: { Attributes: { name } } } = message;

    if (!name.includes('pipeline')) {
      return;
    }

    const coloredName = generateRandomColor(nameColorMap, name);
    console.log('Container with name', coloredName, 'started, attaching log...');

    docker.getContainer(id).attach(
      { stream: true, stdout: true, stderr: true },
      (err, stream) => {
        stream.pipe(pt(`${coloredName} | `)).pipe(process.stdout);
      },
    );
  });

  const stopDieCallback = (message) => {
    const { Actor: { Attributes: { name } } } = message;

    if (!name.includes('pipeline')) {
      return;
    }

    const coloredName = nameColorMap[name] || name;

    console.log('Container with name', coloredName, 'stopped.');
  };

  emitter.on('stop', stopDieCallback);
  emitter.on('die', stopDieCallback);
  emitter.on('destroy', stopDieCallback);
};

const main = async () => {
  // Initialize stack for Lambdas first.

  try {
    await initStack();
  } catch (e) {
    console.error();
    console.error('Could not initialize stack. Are you sure Inframock is running?');
    console.error();
    process.exit(1);
  }

  // Create hook into Docker API.
  const docker = new Docker();
  const emitter = new DockerEvents({
    docker,
  });

  // Create a color map for easy visualization.
  const nameColorMap = {};

  // Attach to existing running containers and new containers, respectively.
  attachToExistingContainers(docker, nameColorMap);
  attachToNewContainers(docker, emitter, nameColorMap);

  console.log('Waiting for Docker events...');
  emitter.start();
};

main();
